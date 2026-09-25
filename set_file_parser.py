import pint

ureg = pint.UnitRegistry(auto_reduce_dimensions=True, system="US")

def read_parameters(filename, section):
    params = {}
    lines = []
    with open(filename) as file:
        lines = file.readlines()

    section_line = None
    for i in range(len(lines)):
        if not lines[i].startswith('#'):
            continue
        section_header = lines[i].lstrip("# ").rstrip()
        if section_header != section:
            continue
        if section_line is not None:
            raise Exception(f"Duplicate Section '{section}' at lines {section_line} and {i}")
        section_line = i
        break
    if section_line is None:
        raise Exception(f"Missing Section '{section}'")

    in_block_header = lines[section_line + 1]
    if not in_block_header.startswith("in:"):
        raise Exception(f"Expected 'in:' block, found '{in_block_header}' at line {section_line + 1}")
    in_block_header = in_block_header[3:]
    in_block_header = in_block_header.lstrip().rstrip()
    if in_block_header != "{":
        raise Exception(f"Expected '{{' to start 'in:' block, found '{in_block_header}' at line {section_line + 1}")

    state = "seek"
    for i in range(section_line + 2, len(lines)):
        line = lines[i].lstrip().rstrip()
        if state == "seek":
            if len(line) == 0:
                continue
            if line == '}':
                break
            parts = line.split('=')
            if len(parts) != 2:
                raise Exception(f"Expected something of the form 'name = value', found '{line}' at line {i}")
            name = parts[0].rstrip()
            value = parts[1].lstrip()
            if name in params:
                raise Exception(f"Duplicate field '{name}' in 'in:' block at line {i}")
            if '"' in value:
                params[name] = value
            else:
                params[name] = ureg(value)
            state = "reason"
        elif state == "reason":
            if not line.startswith("reason:"):
                raise Exception(f"Expected 'reason:' after parameter, found '{line}' at line {i}")
            state = "owner"
        else:
            if not line.startswith("owner:"):
                raise Exception(f"Expected 'owner:' after reason, found '{line}' at line {i}")
            state = "seek"
    return params

def write_outputs(filename, section):
    pass
