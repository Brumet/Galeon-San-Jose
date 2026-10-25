import bpy, sys, os
from mathutils import Vector
src, out, dx, dy_, dz_, sc_ = sys.argv[sys.argv.index("--") + 1:]
bpy.ops.wm.read_factory_settings(use_empty=True); bpy.ops.wm.ply_import(filepath=src); ob = bpy.context.selected_objects[0]
sc = bpy.context.scene; sc.render.engine = 'BLENDER_WORKBENCH'; sh = sc.display.shading
sh.light = 'STUDIO'; sh.color_type = 'SINGLE'; sh.single_color = (0.85, 0.68, 0.32); sh.show_cavity = True; sh.cavity_type = 'BOTH'
sc.world = bpy.data.worlds.new("w")
cam = bpy.data.objects.new("c", bpy.data.cameras.new("c")); sc.collection.objects.link(cam); sc.camera = cam; cam.data.type = 'ORTHO'
bb = [ob.matrix_world @ Vector(c) for c in ob.bound_box]; c = sum(bb, Vector()) / 8
d = Vector((float(dx), float(dy_), float(dz_))).normalized(); cam.location = c + d * 100
cam.rotation_euler = (-d).to_track_quat('-Z', 'Y' if abs(d.z) > 0.9 else 'Z').to_euler(); cam.data.ortho_scale = float(sc_)
sc.render.resolution_x, sc.render.resolution_y = 1100, 900; sc.render.filepath = out; bpy.ops.render.render(write_still=True)
