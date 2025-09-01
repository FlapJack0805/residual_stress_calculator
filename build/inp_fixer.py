import numpy as np, re

inp = open("calculix.inp").read()

# Parse nodes
nodes = {}
reading_nodes = False
for line in inp.splitlines():
    if line.upper().startswith("*NODE"):
        reading_nodes = True; continue
    if reading_nodes:
        if line.strip().startswith("*"): break
        if not line.strip(): continue
        parts = [p.strip() for p in line.split(",")]
        if len(parts) >= 4 and parts[0].replace("-","").isdigit():
            nid = int(parts[0]); x,y,z = map(float, parts[1:4])
            nodes[nid] = np.array([x,y,z])

# Parse C3D4 elements
elems = []
lines = inp.splitlines()
i = 0
while i < len(lines):
    line = lines[i]
    if line.upper().startswith("*ELEMENT") and "C3D4" in line.upper():
        i += 1
        while i < len(lines) and not lines[i].strip().startswith("*"):
            if lines[i].strip():
                parts = [p.strip() for p in lines[i].split(",")]
                if len(parts) >= 5:
                    eid = int(parts[0])
                    n = list(map(int, parts[1:5]))
                    elems.append((eid,n))
            i += 1
        continue
    i += 1

def svol(n):
    p1,p2,p3,p4 = [nodes[j] for j in n]
    M = np.column_stack((p2-p1, p3-p1, p4-p1))
    return np.linalg.det(M)/6.0

neg, tiny, dup = [], [], []
for eid, n in elems:
    if len(set(n)) < 4:
        dup.append(eid); continue
    V = svol(n)
    # size-aware sliver test
    ps = [nodes[j] for j in n]
    L = max(np.linalg.norm(ps[a]-ps[b]) for a in range(4) for b in range(a+1,4))
    if V < 0: neg.append(eid)
    elif abs(V) < 1e-6*(L**3): tiny.append(eid)

print("neg =", len(neg), "tiny =", len(tiny), "dup =", len(dup))

# Write a patched file flipping orientation of negative tets (swap first two nodes)
out = []
in_c3d4 = False
negset = set(neg)
for line in lines:
    if line.upper().startswith("*ELEMENT") and "C3D4" in line.upper():
        in_c3d4 = True; out.append(line); continue
    if in_c3d4 and line.strip().startswith("*"):
        in_c3d4 = False; out.append(line); continue
    if in_c3d4 and line.strip():
        parts = [p.strip() for p in line.split(",")]
        if parts[0].lstrip("-").isdigit():
            eid = int(parts[0])
            if eid in negset:
                parts[1], parts[2] = parts[2], parts[1]  # flip
                out.append(",".join(parts))
            else:
                out.append(line)
        else:
            out.append(line)
    else:
        out.append(line)

open("calculix.inp","w").write("\n".join(out))
