# abaqus_tester.py
# Usage:
#   abaqus cae noGUI=abaqus_tester.py -- init --case case_000 --seed 7
#   abaqus cae noGUI=abaqus_tester.py -- cut  --case case_000 --cut_labels cut_001.txt --job CUT_001
#
# Notes:
# - Cuts are done by ELEMENT deactivation (Abaqus doesn't "deactivate nodes" for analysis).
# - Each cut is a separate job restarted from the previous job's end state.
# - Exports:
#     case_000/BASE.odb, disp_BASE.csv, mesh_BASE.vtk, state.json
#     case_000/step_###/CUT_###.odb, disp_CUT_###.csv, mesh_CUT_###.vtk, DONE.txt

from abaqus import *
from abaqusConstants import *
import regionToolset
import mesh as abqmesh
import os, sys, json, csv, math, random

# -----------------------
# FS helpers
# -----------------------

def ensure_dir(p):
    if p and not os.path.exists(p):
        os.makedirs(p)

def abs_path(p):
    return os.path.abspath(p)

# -----------------------
# Simple CLI parse
# -----------------------

def parse_args(argv):
    # argv already excludes the Abaqus launcher args; includes subcommand at argv[0]
    if not argv:
        raise RuntimeError("Need subcommand: init | cut")

    cmd = argv[0]
    args = argv[1:]

    # common
    out = {'cmd': cmd}

    i = 0
    while i < len(args):
        a = args[i]
        if a == '--case':
            out['case'] = args[i+1]; i += 2
        elif a == '--seed':
            out['seed'] = int(args[i+1]); i += 2
        elif a == '--job':
            out['job'] = args[i+1]; i += 2
        elif a == '--cut_labels':
            out['cut_labels'] = args[i+1]; i += 2

        # geometry knobs for init (optional)
        elif a == '--L':
            out['L'] = float(args[i+1]); i += 2
        elif a == '--W':
            out['W'] = float(args[i+1]); i += 2
        elif a == '--H':
            out['H'] = float(args[i+1]); i += 2
        elif a == '--elem':
            out['elem'] = float(args[i+1]); i += 2
        elif a == '--pockets':
            out['pockets'] = int(args[i+1]); i += 2
        elif a == '--bumps':
            out['bumps'] = int(args[i+1]); i += 2

        else:
            raise RuntimeError("Unknown arg: %s" % a)

    return out

# -----------------------
# Model building for BASE
# -----------------------

def bbox_of_part(part):
    bb = part.getBoundingBox()
    lo = bb['low']; hi = bb['high']
    return (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2])

def clamp(x, a, b):
    return max(a, min(b, x))

def build_block_with_pockets(model, part_name, L, W, H, num_pockets, seed):
    random.seed(seed)

    s = model.ConstrainedSketch(name='__profile__', sheetSize=10*max(L, W))
    s.rectangle(point1=(0.0, 0.0), point2=(L, W))

    part = model.Part(name=part_name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    part.BaseSolidExtrude(sketch=s, depth=H)
    del s

    margin = 0.08 * min(L, W)
    for _ in range(num_pockets):
        px = random.uniform(margin, L - margin)
        py = random.uniform(margin, W - margin)
        pw = random.uniform(0.12*L, 0.35*L)
        ph = random.uniform(0.12*W, 0.35*W)

        x1 = clamp(px - pw/2.0, margin, L - margin)
        x2 = clamp(px + pw/2.0, margin, L - margin)
        y1 = clamp(py - ph/2.0, margin, W - margin)
        y2 = clamp(py + ph/2.0, margin, W - margin)

        depth = random.uniform(0.15*H, 0.75*H)

        sk = model.ConstrainedSketch(name='__pocket__', sheetSize=10*max(L, W))
        sk.rectangle(point1=(x1, y1), point2=(x2, y2))

        top_face = part.faces.findAt(((L*0.5, W*0.5, H),))[0]
        up_edge  = part.edges.findAt(((L, W*0.5, H),))[0]

        part.CutExtrude(sketchPlane=top_face,
                        sketchPlaneSide=SIDE1,
                        sketchUpEdge=up_edge,
                        sketch=sk,
                        flipExtrudeDirection=ON,
                        depth=depth)
        del sk

    return part

def assign_material_section(model, part, E=70e3, nu=0.33, alpha=23e-6):
    mat = model.Material(name='Mat')
    mat.Elastic(table=((E, nu),))
    mat.Expansion(table=((alpha,),))
    sec = model.HomogeneousSolidSection(name='Sec', material='Mat', thickness=None)
    region = regionToolset.Region(cells=part.cells)
    part.SectionAssignment(region=region, sectionName='Sec')
    return sec

def mesh_part(part, elem_size):
    part.seedPart(size=elem_size, deviationFactor=0.1, minSizeFactor=0.1)
    et = abqmesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)
    part.setElementType(regions=(part.cells,), elemTypes=(et,))
    part.generateMesh()

def create_instance(model, part):
    asm = model.rootAssembly
    asm.DatumCsysByDefault(CARTESIAN)
    inst = asm.Instance(name='Part-1', part=part, dependent=ON)
    return inst

def apply_min_rbm_constraints(model, inst):
    asm = model.rootAssembly
    part = inst.part
    xmin,xmax,ymin,ymax,zmin,zmax = bbox_of_part(part)

    p1 = (xmin, ymin, zmin)
    p2 = (xmax, ymin, zmin)
    p3 = (xmin, ymax, zmin)

    v1 = inst.vertices.findAt((p1,))
    v2 = inst.vertices.findAt((p2,))
    v3 = inst.vertices.findAt((p3,))

    asm.Set(name='RBM_1', vertices=v1)
    asm.Set(name='RBM_2', vertices=v2)
    asm.Set(name='RBM_3', vertices=v3)

    model.DisplacementBC('BC_RBM_1', 'Initial', asm.sets['RBM_1'], u1=0.0, u2=0.0, u3=0.0)
    model.DisplacementBC('BC_RBM_2', 'Initial', asm.sets['RBM_2'], u2=0.0, u3=0.0)
    model.DisplacementBC('BC_RBM_3', 'Initial', asm.sets['RBM_3'], u3=0.0)

def gaussian_temp(x,y,z, bumps):
    T = 0.0
    for (cx,cy,cz,A,sig) in bumps:
        dx=x-cx; dy=y-cy; dz=z-cz
        r2 = dx*dx + dy*dy + dz*dz
        T += A * math.exp(-r2/(2.0*sig*sig))
    return T

def make_random_bumps(part, n_bumps, seed, A_range=(40.0, 180.0), sigma_frac=(0.12, 0.35)):
    random.seed(seed)
    xmin,xmax,ymin,ymax,zmin,zmax = bbox_of_part(part)
    L = xmax-xmin; W=ymax-ymin; H=zmax-zmin
    m = min(L,W,H)
    bumps=[]
    for _ in range(n_bumps):
        cx = random.uniform(xmin+0.1*L, xmax-0.1*L)
        cy = random.uniform(ymin+0.1*W, ymax-0.1*W)
        cz = random.uniform(zmin+0.1*H, zmax-0.1*H)
        A  = random.uniform(A_range[0], A_range[1])
        sig = random.uniform(sigma_frac[0]*m, sigma_frac[1]*m)
        bumps.append((cx,cy,cz,A,sig))
    return bumps

def apply_temperature_field(model, inst, bumps, seed_tag):
    asm = model.rootAssembly
    asm.Set(name='ALL_NODES', nodes=inst.nodes)

    data=[]
    for nd in inst.nodes:
        x,y,z = nd.coordinates
        T = gaussian_temp(x,y,z,bumps)
        data.append((nd.label, T))

    field_name = 'Tfield_%s' % seed_tag
    model.DiscreteField(name=field_name, defaultValues=(0.0,), data=data)

    model.Temperature(name='PredefTemp',
                      createStepName='ProcessResidual',
                      region=asm.sets['ALL_NODES'],
                      distributionType=FIELD,
                      field=field_name)

def zero_temperature_in_relax(model):
    asm = model.rootAssembly
    model.Temperature(name='TempToZero',
                      createStepName='RelaxResidual',
                      region=asm.sets['ALL_NODES'],
                      distributionType=UNIFORM,
                      magnitudes=(0.0,))

def add_steps(model):
    model.StaticStep(name='ProcessResidual', previous='Initial', nlgeom=OFF)
    model.StaticStep(name='RelaxResidual', previous='ProcessResidual', nlgeom=OFF)

def request_outputs(model):
    # Make sure we have U in ODB; STATUS will not exist in BASE (no removals yet)
    model.fieldOutputRequests['F-Output-1'].setValues(variables=('S','U','E'))

def inject_restart_write_keyword(model, step_name='RelaxResidual'):
    """
    Ensure restart data is written at end of RelaxResidual so future jobs can *RESTART, READ.
    Most reliable way across versions: inject keyword before the *END STEP of that step.
    """
    model.keywordBlock.synchVersions()

    lines = list(model.keywordBlock.sieBlocks)

    # Find the RelaxResidual step block and insert *RESTART, WRITE before its *END STEP
    # This is simplistic but usually works because CAE emits:
    # *Step, name=RelaxResidual ...
    # ...
    # *End Step
    i_step = None
    i_end  = None
    for i,ln in enumerate(lines):
        l = ln.strip().lower()
        if l.startswith('*step') and ('relaxresidual' in l):
            i_step = i
        if i_step is not None and l.startswith('*end step'):
            i_end = i
            break

    if i_step is None or i_end is None:
        # fallback: append at end (better than nothing)
        lines.append('** Restart write (fallback)\n')
        lines.append('*RESTART, WRITE, FREQUENCY=1\n')
    else:
        ins = i_end
        lines.insert(ins, '*RESTART, WRITE, FREQUENCY=1\n')

    model.keywordBlock.sieBlocks = tuple(lines)

# -----------------------
# ODB exporters (truth)
# -----------------------

def export_disp_csv(odb_path, step_name, out_csv):
    from odbAccess import openOdb
    odb=openOdb(odb_path, readOnly=True)
    step=odb.steps[step_name]
    frame=step.frames[-1]
    U=frame.fieldOutputs['U']
    inst_name = odb.rootAssembly.instances.keys()[0]
    inst = odb.rootAssembly.instances[inst_name]
    node_coords = {n.label: n.coordinates for n in inst.nodes}

    with open(out_csv,'w') as f:
        w=csv.writer(f)
        w.writerow(['nodeLabel','x','y','z','u1','u2','u3'])
        for v in U.values:
            lab=v.nodeLabel
            x,y,z=node_coords[lab]
            u1,u2,u3=v.data
            w.writerow([lab,x,y,z,u1,u2,u3])
    odb.close()

def export_active_mesh_vtk(odb_path, step_name, out_vtk):
    """
    Exports volume mesh of active elements. If STATUS exists (cuts), use it.
    If not (BASE), treat all as active.
    """
    from odbAccess import openOdb
    odb=openOdb(odb_path, readOnly=True)
    step=odb.steps[step_name]
    frame=step.frames[-1]
    inst_name = odb.rootAssembly.instances.keys()[0]
    inst = odb.rootAssembly.instances[inst_name]

    status_map={}
    if 'STATUS' in frame.fieldOutputs.keys():
        for v in frame.fieldOutputs['STATUS'].values:
            status_map[v.elementLabel]=v.data

    active=[]
    for e in inst.elements:
        if status_map.get(e.label, 1.0) > 0.5:
            active.append(e)

    used=set()
    for e in active:
        for nl in e.connectivity:
            used.add(nl)

    node_coord={n.label:n.coordinates for n in inst.nodes if n.label in used}
    labs=sorted(node_coord.keys())
    vtk_id={lab:i for i,lab in enumerate(labs)}

    VTK_QUADRATIC_TET=24  # 10-node tet

    with open(out_vtk,'w') as f:
        f.write('# vtk DataFile Version 3.0\n')
        f.write('Active mesh %s\n' % step_name)
        f.write('ASCII\n')
        f.write('DATASET UNSTRUCTURED_GRID\n')
        f.write('POINTS %d float\n' % len(labs))
        for lab in labs:
            x,y,z=node_coord[lab]
            f.write('%g %g %g\n' % (x,y,z))

        total_ints = len(active) * (1 + 10)
        f.write('CELLS %d %d\n' % (len(active), total_ints))
        for e in active:
            conn=[vtk_id[nl] for nl in e.connectivity]
            f.write('10 ' + ' '.join(str(i) for i in conn) + '\n')

        f.write('CELL_TYPES %d\n' % len(active))
        for _ in active:
            f.write('%d\n' % VTK_QUADRATIC_TET)

    odb.close()

# -----------------------
# Restart-cut input writer
# -----------------------

def read_elem_labels(path):
    labels=[]
    with open(path,'r') as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith('#'): continue
            for tok in line.replace(',',' ').split():
                labels.append(int(tok))
    if not labels:
        raise RuntimeError('No element labels in %s' % path)
    return labels

def write_restart_cut_inp(inp_path, prev_job, prev_step, labels):
    with open(inp_path, 'w') as f:
        f.write('*HEADING\n')
        f.write('** Restart cut\n')
        f.write('*RESTART, READ, JOB=%s, STEP=%s, INC=0\n' % (prev_job, prev_step))
        f.write('*ELSET, ELSET=TO_REMOVE\n')
        for i,lab in enumerate(labels):
            f.write('%d' % lab)
            if i == len(labels)-1:
                f.write('\n')
            elif (i+1) % 16 == 0:
                f.write(',\n')
            else:
                f.write(', ')
        f.write('*STEP\n')
        f.write('*STATIC\n')
        f.write('*MODEL CHANGE, REMOVE\n')
        f.write('TO_REMOVE\n')
        f.write('*OUTPUT, FIELD, FREQ=1\n')
        f.write('*NODE OUTPUT\n')
        f.write('U\n')
        f.write('*ELEMENT OUTPUT\n')
        f.write('STATUS\n')
        f.write('*END STEP\n')

# -----------------------
# Actions: init / cut
# -----------------------

def action_init(case_dir, seed, L, W, H, elem_size, pockets, bumps, job='BASE'):
    ensure_dir(case_dir)
    case_dir = abs_path(case_dir)

    # Run BASE in the case directory so BASE.res, BASE.odb, etc live there.
    os.chdir(case_dir)

    # Fresh model
    mname = 'MODEL_BASE'
    if mname in mdb.models:
        del mdb.models[mname]
    model = mdb.Model(name=mname)

    part = build_block_with_pockets(model, 'Part', L, W, H, pockets, seed)
    assign_material_section(model, part)

    mesh_part(part, elem_size)
    inst = create_instance(model, part)

    add_steps(model)
    apply_min_rbm_constraints(model, inst)

    b = make_random_bumps(part, n_bumps=bumps, seed=seed+1000)
    apply_temperature_field(model, inst, b, seed_tag=str(seed))
    zero_temperature_in_relax(model)

    request_outputs(model)

    # Make sure restart is written at end of RelaxResidual
    inject_restart_write_keyword(model, step_name='RelaxResidual')

    # Submit BASE job
    j = mdb.Job(name=job, model=mname, type=ANALYSIS,
                numCpus=4, numDomains=4,
                description='BASE residual state seed=%d' % seed)
    j.submit(consistencyChecking=OFF)
    j.waitForCompletion()

    odb = job + '.odb'
    export_disp_csv(odb, 'RelaxResidual', 'disp_%s.csv' % job)
    export_active_mesh_vtk(odb, 'RelaxResidual', 'mesh_%s.vtk' % job)

    # Save state.json for restart loop
    state = {
        "case_dir": case_dir,
        "current_job": job,
        "current_step": "RelaxResidual",
        "cut_index": 0
    }
    with open(os.path.join(case_dir, 'state.json'), 'w') as f:
        json.dump(state, f, indent=2)

    with open(os.path.join(case_dir, 'DONE_INIT.txt'), 'w') as f:
        f.write('ok\n')

def action_cut(case_dir, cut_labels_file, job):
    case_dir = abs_path(case_dir)
    state_path = os.path.join(case_dir, 'state.json')
    if not os.path.exists(state_path):
        raise RuntimeError("Missing state.json in %s (run init first)" % case_dir)

    with open(state_path,'r') as f:
        state = json.load(f)

    prev_job  = state['current_job']
    prev_step = state['current_step']
    cut_index = int(state.get('cut_index', 0)) + 1

    labels = read_elem_labels(cut_labels_file)

    step_dir = os.path.join(case_dir, 'step_%03d' % cut_index)
    ensure_dir(step_dir)

    # Keep a copy of labels used for this step
    local_labels = os.path.join(step_dir, 'cut_elems.txt')
    with open(cut_labels_file,'r') as src, open(local_labels,'w') as dst:
        dst.write(src.read())

    os.chdir(step_dir)

    inp = job + '.inp'
    write_restart_cut_inp(inp, prev_job=prev_job, prev_step=prev_step, labels=labels)

    # IMPORTANT: run Abaqus Standard from within step_dir, but restart reads BASE.res etc from case_dir.
    # Abaqus finds restart files by JOB name; so keep case_dir on the path by running in case_dir?
    #
    # Easiest reliable method: run the job from case_dir, but write outputs into step_dir by changing workdir.
    # Abaqus doesn't have a great "workdir" flag everywhere, so we do this:
    # - symlink/copy restart files into step_dir OR
    # - run from case_dir and use relative output paths
    #
    # We'll do the simplest robust thing: run from case_dir, but keep the cut artifacts in step_dir by prefixing job name.
    #
    # Instead, we assume BASE.* files are in case_dir and step_dir is inside case_dir.
    # Restart read should work if you run from inside step_dir (Abaqus will look in current dir).
    # So we copy required restart files from case_dir into step_dir (small overhead, robust).

    # Copy all files that start with prev_job (BASE.*) into step_dir
    for fn in os.listdir(case_dir):
        if fn.startswith(prev_job + '.'):
            src = os.path.join(case_dir, fn)
            dst = os.path.join(step_dir, fn)
            if not os.path.exists(dst):
                try:
                    with open(src,'rb') as a, open(dst,'wb') as b:
                        b.write(a.read())
                except:
                    pass

    cmd = 'abaqus job=%s input=%s interactive' % (job, inp)
    rc = os.system(cmd)
    if rc != 0:
        raise RuntimeError("Cut job failed (rc=%d). See %s/%s.msg" % (rc, step_dir, job))

    odb = job + '.odb'
    export_disp_csv(odb, 'Step-1', 'disp_%s.csv' % job)
    export_active_mesh_vtk(odb, 'Step-1', 'mesh_%s.vtk' % job)

    with open(os.path.join(step_dir, 'DONE.txt'), 'w') as f:
        f.write('ok\n')

    # Update global state for next cut
    state['current_job']  = job
    state['current_step'] = 'Step-1'
    state['cut_index']    = cut_index
    with open(state_path,'w') as f:
        json.dump(state, f, indent=2)

# -----------------------
# Entry
# -----------------------

def main():
    argv = sys.argv
    argv = argv[argv.index('--')+1:] if '--' in argv else []
    a = parse_args(argv)

    cmd = a.get('cmd')
    case_dir = a.get('case', 'case_000')

    if cmd == 'init':
        seed = a.get('seed', 1)
        job  = a.get('job', 'BASE')
        # defaults (override via flags if you want)
        L = a.get('L', 120.0)
        W = a.get('W', 70.0)
        H = a.get('H', 35.0)
        elem_size = a.get('elem', 5.0)
        pockets = a.get('pockets', 4)
        bumps = a.get('bumps', 5)
        action_init(case_dir, seed, L, W, H, elem_size, pockets, bumps, job=job)

    elif cmd == 'cut':
        job = a.get('job', None)
        cut_labels = a.get('cut_labels', None)
        if not job or not cut_labels:
            raise RuntimeError("cut needs --job and --cut_labels")
        action_cut(case_dir, cut_labels, job)

    else:
        raise RuntimeError("Unknown command: %s" % cmd)

if __name__ == '__main__':
    main()
