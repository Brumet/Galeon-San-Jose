import bpy, sys, os, math
from mathutils import Vector, Euler
src, out, tag = sys.argv[sys.argv.index("--") + 1:]
bpy.ops.wm.read_factory_settings(use_empty=True)
(bpy.ops.wm.ply_import if src.endswith(".ply") else bpy.ops.wm.stl_import)(filepath=src)
sc = bpy.context.scene
sc.render.engine = 'BLENDER_WORKBENCH'
sh = sc.display.shading
sh.light = 'STUDIO'; sh.color_type = 'SINGLE'; sh.single_color = (0.85, 0.68, 0.32)
sh.show_cavity = True; sh.cavity_type = 'BOTH'
sc.world = bpy.data.worlds.new("w")
cam = bpy.data.objects.new("cam", bpy.data.cameras.new("cam")); sc.collection.objects.link(cam); sc.camera = cam
cam.data.type = 'ORTHO'
c = Vector((110, 0, 21.8))
views = {
 "abajo": (Vector((0,0,-300)), Euler((math.pi,0,0)), 230, (1800,500)),
 "lado_menosY": (Vector((0,-300,0)), Euler((math.pi/2,0,0)), 230, (1800,500)),
 "lado_masY": (Vector((0,300,0)), Euler((math.pi/2,0,math.pi)), 230, (1800,500)),
 "tresc_abajo": (Vector((0,-212,-212)), Euler((math.radians(135),0,0)), 230, (1800,600)),
 "munon_masY": (Vector((8,100,-105)), Euler((math.radians(-135),0,0)), 40, (1200,1000)),
 "munon_menosY": (Vector((8,-100,-105)), Euler((math.radians(135),0,0)), 40, (1200,1000)),
 "culata_abajo": (Vector((80,-150,-150)), Euler((math.radians(135),0,0)), 60, (1400,900)),
}
for n, (off, rot, s, res) in views.items():
    cam.location = c + off; cam.rotation_euler = rot; cam.data.ortho_scale = s
    sc.render.resolution_x, sc.render.resolution_y = res
    sc.render.filepath = os.path.join(out, f"{tag}_{n}.png")
    bpy.ops.render.render(write_still=True); print("rendered", n, flush=True)
