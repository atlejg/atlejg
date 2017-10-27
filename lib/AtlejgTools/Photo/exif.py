from PIL import Image
from PIL.ExifTags import TAGS

def get_exif(fname):
   ret = {}
   i = Image.open(fname)
   info = i._getexif()
   for tag, value in info.items():
      decoded = TAGS.get(tag, tag)
      ret[decoded] = value
   return ret

fname = 'TONE.jpg'
exif = get_exif(fname)
for (key, value) in exif.items():
   # dont print long lists..
   if type(value) is int or len(value) < 30:
      print key, value

print 'aperture is:', exif['ApertureValue'][0] / float(exif['ApertureValue'][1])

