import bpy, sys, os, math
from mathutils import Vector
src, out, tag = sys.argv[sys.argv.index("--") + 1:]
bpy.ops.wm.read_factory_settings(use_empty=True)
bpy.ops.wm.ply_import(filepath=src)
ob = bpy.context.selected_objects[0]
sc = bpy.context.scene; sc.render.engine = 'BLENDER_WORKBENCH'
sh = sc.display.shading; sh.light = 'STUDIO'; sh.show_cavity = True; sh.cavity_type = 'BOTH'
try: sh.color_type = 'VERTEX'
except Exception: sh.color_type = 'ATTRIBUTE'
sc.world = bpy.data.worlds.new("w")
cam = bpy.data.objects.new("cam", bpy.data.cameras.new("cam")); sc.collection.objects.link(cam); sc.camera = cam
cam.data.type = 'ORTHO'
bb = [ob.matrix_world @ Vector(c) for c in ob.bound_box]
c = sum(bb, Vector()) / 8; size = max((max(v[i] for v in bb) - min(v[i] for v in bb)) for i in range(3))
dirs = {"frente": (0,-1,0), "atras": (0,1,0), "arriba": (0,0,1), "izq": (-1,0,0), "der": (1,0,0), "abajo": (0,0,-1), "diag1": (-1,-1,0.8), "diag2": (1,1,0.8)}
sc.render.resolution_x, sc.render.resolution_y = 1100, 900
for n, d in dirs.items():
    d = Vector(d).normalized(); cam.location = c + d * 200
    cam.rotation_euler = (-d).to_track_quat('-Z', 'Z' if abs(d.z) < 0.99 else 'Y').to_euler()
    cam.data.ortho_scale = size * 1.1
    sc.render.filepath = os.path.join(out, f"{tag}_{n}.png"); bpy.ops.render.render(write_still=True)
