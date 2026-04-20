import pyvista as pv
import sys

if len(sys.argv) != 3:
    print("Usage: python convert_stl_to_vtk.py input.stl output.vtk")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

mesh = pv.read(input_file)

print("Bounds:", mesh.bounds)
print("Number of faces:", mesh.n_cells)

mesh.save(output_file)

print(f"Saved to {output_file}")
