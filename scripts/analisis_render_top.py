import bpy, sys, os
from mathutils import Vector
src, out = sys.argv[sys.argv.index("--") + 1:]
bpy.ops.wm.read_factory_settings(use_empty=True); bpy.ops.wm.ply_import(filepath=src); ob = bpy.context.selected_objects[0]
sc = bpy.context.scene; sc.render.engine = 'BLENDER_WORKBENCH'; sh = sc.display.shading
sh.light = 'STUDIO'; sh.color_type = 'SINGLE'; sh.single_color = (0.85, 0.68, 0.32); sh.show_cavity = True; sh.cavity_type = 'BOTH'
sc.world = bpy.data.worlds.new("w")
cam = bpy.data.objects.new("c", bpy.data.cameras.new("c")); sc.collection.objects.link(cam); sc.camera = cam; cam.data.type = 'ORTHO'
bb = [ob.matrix_world @ Vector(c) for c in ob.bound_box]; c = sum(bb, Vector()) / 8
cam.location = c + Vector((0, 0, 100)); cam.rotation_euler = (0, 0, 0); cam.data.ortho_scale = 17
sc.render.resolution_x, sc.render.resolution_y = 1100, 1000; sc.render.filepath = out; bpy.ops.render.render(write_still=True)
