#

from pdf2image import convert_from_path

import sys

in_file = sys.argv[1]
out_file = in_file.split(".")[0] + ".jpg"

pages = convert_from_path(in_file, 200)

for page in pages:
    page.save(out_file, 'JPEG') 
