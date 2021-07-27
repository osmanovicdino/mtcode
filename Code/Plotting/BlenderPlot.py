import bpy
import bmesh
from itertools import product

def create_base_sphere(context, name, size):
    # Create an object
    mesh = bpy.data.meshes.new(name)
    obj = bpy.data.objects.new(name, mesh)

    # Add a cube
    bm = bmesh.new()
    bm.from_mesh(mesh)
    # bmesh.ops.create_cube(bm, size=size, calc_uvs=False)
    bpy.ops.mesh.primitive_uv_sphere_add(segments=16, ring_count=8, location=(float(0), float(0), float(0)),  radius=0.5)
    bm.to_mesh(mesh)
    bm.free()

    # Link the object to the scene and set it active
    context.scene.collection.objects.link(obj)
    context.view_layer.objects.active = obj

    return obj


def make_result_object(context, name, cube, vertices):
    mesh = bpy.data.meshes.new(name)
    obj = bpy.data.objects.new(name, mesh)

    # Associate the vertices
    obj.data.from_pydata(vertices, [], [])

    # Link the object to the scene and set it active
    context.scene.collection.objects.link(obj)
    context.view_layer.objects.active = obj

    # Make the object parent of the cube
    cube.parent = obj
    # Make the object dupliverts
    obj.instance_type = 'VERTS'

    return obj
 

#  import bpy

# text_file = open("/home/dino/External/PhaseDiagramDesign6/den=0.005_i1=12._i2=12._i3=4._i4=60._m1=16000_m2=16000_rate=0.001/posden=0.005_int1=12_int2=12_int3=4_int4=60_br=0.001num_anti=16000num_inv=16000_i=00000.csv", "r")
# lines = []

# #Read in contents to a list
# for line in text_file:
#     lines.append(line.strip())

# it = 0
# #create spheres
# for e in lines:
#     print(it)
#     it = it + 1
#     temp = e.split(',')
#     bpy.ops.mesh.primitive_uv_sphere_add(segments=16, ring_count=8 ,location=(float(temp[0]), float(temp[1]), float(temp[2])),  radius=0.5)
#     bpy.ops.object.shade_smooth()

# text_file.close()


# import bpy

# for line in open("spheres.txt"):
#     // you have to parse "x", "y", "z" and "r" from the variable "line"
#     x, y, z, r = line.split(",")
#     // making the sphere
#     bpy.ops.mesh.primitive_ico_sphere_add(size=r, location=(x,y,z)) 
