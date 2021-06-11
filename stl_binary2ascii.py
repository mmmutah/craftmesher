import struct
import argparse

def unpack (f, sig, l, fb):
    s = f.read (l)
    fb.append(s)
    return struct.unpack(sig, s)

def read_triangle(f,normals,points,triangles,bytecount,fb):
    n = unpack(f,"<3f", 12, fb)
    p1 = unpack(f,"<3f", 12, fb)
    p2 = unpack(f,"<3f", 12, fb)
    p3 = unpack(f,"<3f", 12, fb)
    b = unpack(f,"<h", 2, fb)

    normals.append(n)
    l = len(points)
    points.append(p1)
    points.append(p2)
    points.append(p3)
    triangles.append((l, l+1, l+2))
    bytecount.append(b[0])

def read_length(f):
    length = struct.unpack("@i", f.read(4))
    return length[0]

def read_header(f):
    f.seek(f.tell()+80)

def write_as_ascii(outfilename,normals,points,triangles,bytecount,fb):

    f = open(outfilename, "w")
    f.write ("solid\n")
    for n  in range(len(triangles)):
        f.write (" facet normal {} {} {}\n".format(normals[n][0],normals[n][1],normals[n][2]))
        f.write ("  outer loop\n")
        f.write ("   vertex {} {} {}\n".format(points[triangles[n][0]][0],points[triangles[n][0]][1],points[triangles[n][0]][2]))
        f.write ("   vertex {} {} {}\n".format(points[triangles[n][1]][0],points[triangles[n][1]][1],points[triangles[n][1]][2]))
        f.write ("   vertex {} {} {}\n".format(points[triangles[n][2]][0],points[triangles[n][2]][1],points[triangles[n][2]][2]))
        f.write ("  endloop\n")
        f.write (" endfacet\n")
    f.write ("endsolid\n")
    f.close()

def convert(filename):

    normals = []
    points = []
    triangles = []
    bytecount = []

    fb = [] # debug list
    
    infilename = filename
    outfilename = 'ascii' + filename

    try:
        f = open ( infilename, "rb")

        read_header(f)
        l = read_length(f)
        try:
            while True:
                read_triangle(f,normals,points,triangles,bytecount,fb)
        except Exception as e:
            mm = 0

        write_as_ascii(outfilename,normals,points,triangles,bytecount,fb)

    except Exception as e:
        print(e)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process the job name.')
    parser.add_argument('-j', action="store")
    args, unknown = parser.parse_known_args()

    jobname = str(args.j)
    print(jobname)
    convert(jobname)
