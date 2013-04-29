import bpy

file = '/Users/mattvisco/Documents/School/CS184/fluid_simulation/output.txt'
f = open(file)
frame_count = 0

bpy.ops.object.select_all()
bpy.data.objects["Cube"].select = True
bpy.ops.object.delete()

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
        bpy.data.objects["Camera"].select = False
        bpy.data.objects["Lampe"].select = False
        bpy.ops.object.delete() # Deletes selected objects
        bpy.data.objects["Camera"].select = True
        bpy.data.objects["Lamp"].select = True
    else:
        bpy.ops.object.metaball_add(type='BALL', location=(float(pos[0]),float(pos[1]),float(pos[2])))


bpy.data.scenes[0].render.fps = 24
bpy.data.scenes[0].render.resolution_x = 320
bpy.data.scenes[0].render.resolution_y = 240
bpy.data.scenes[0].frame_end = frame_count
bpy.data.scenes[0].render.filepath = '/Users/mattvisco/Documents/School/CS184/fluid_simulation/simulation'

# This is for image output
bpy.data.scenes[0].render.image_settings.file_format = 'JPEG'
bpy.ops.render.render( write_still=True )

# This is for movie output
#bpy.data.scenes[0].render.is_movie_format = True #Maybe?
#bpy.ops.render.render(animation = True)
