
def load_header(filename):
    file = open(filename)
    lines = []
    for line in file.readlines():
        if line[0:2] == '#@':
            hline = line[2:].lstrip()
            lines.append(hline)

    g = {}
    l = {}
    for line in lines:
        exec(line, g, l)

    return l
