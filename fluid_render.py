import bpy

file = 'output.txt'
f = open(file)
frame_count = 0

bpy.data.objects["Camera"].select = True
bpy.ops.transform.translate(value = (0,0,-30))
bpy.ops.transform.rotate(value = (25), axis = (1,0,0))
bpy.ops.transform.rotate(value = (-35), axis = (0,1,0))
bpy.data.objects["Camera"].select = False

for line in f:
    pos = line.split()
    if (line == "FRAME" or line == "FRAME\n"):
        frame_count += 1
        bpy.ops.object.select_all() # DESELECT
        bpy.ops.object.select_all() # SELECT
        bpy.context.scene.frame_current = frame_count
        bpy.ops.anim.keyframe_insert(type = "LocRotScale") # look into other types
        bpy.ops.object.delete() # DELETE
    else:
        bpy.ops.object.metaball_add(type='BALL', location=(pos[0],pos[1],pos[2]))


bpy.data.scenes[0].render.fps = 24
bpy.data.scenes[0].render.resolution_x = 320
bpy.data.scenes[0].render.resolution_y = 240
bpy.data.scenes[0].frame_end = frame_count
bpy.data.scenes[0].render.file_format = "AVI_RAW"
bpy.ops.render.render(animation = True)
