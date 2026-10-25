import bpy, sys, os, math
from mathutils import Vector
src, out, tag, side = sys.argv[sys.argv.index("--") + 1:]
bpy.ops.wm.read_factory_settings(use_empty=True); bpy.ops.wm.ply_import(filepath=src)
ob = bpy.context.selected_objects[0]
sc = bpy.context.scene; sc.render.engine = 'BLENDER_WORKBENCH'
sh = sc.display.shading; sh.light = 'STUDIO'; sh.color_type = 'SINGLE'; sh.single_color = (0.85, 0.68, 0.32); sh.show_cavity = True; sh.cavity_type = 'BOTH'
sc.world = bpy.data.worlds.new("w")
cam = bpy.data.objects.new("cam", bpy.data.cameras.new("cam")); sc.collection.objects.link(cam); sc.camera = cam; cam.data.type = 'ORTHO'
bb = [ob.matrix_world @ Vector(c) for c in ob.bound_box]; c = sum(bb, Vector()) / 8
sg = 1 if side == "masY" else -1
for n, d in {"lado": (0, sg, 0.15), "cabeza": (-1, sg * 0.5, 0.3), "arriba": (0.0, sg * 0.3, 1)}.items():
    d = Vector(d).normalized(); cam.location = c + d * 100
    cam.rotation_euler = (-d).to_track_quat('-Z', 'Z' if abs(d.z) < 0.95 else 'Y').to_euler(); cam.data.ortho_scale = 22 if n != "cabeza" else 14
    sc.render.resolution_x, sc.render.resolution_y = 1000, 800
    sc.render.filepath = os.path.join(out, f"{tag}_{n}.png"); bpy.ops.render.render(write_still=True)
