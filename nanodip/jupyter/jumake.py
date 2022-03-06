#!/usr/bin/env python3

"""
Run this script to parses the whole nanodip code and merge into one single
jupyter notebook file whithin this directory.
"""

import nbformat as nbf
import re
from distutils.dir_util import copy_tree

START_EXTERNAL = '# start_external_modules\n'
END_EXTERNAL = '# end_external_modules\n'
START_INTERNAL = '# start_internal_modules\n'
END_INTERNAL = '# end_internal_modules\n'

def external_imports(modules):
    """Return all external module imports used."""
    def external_modules(content):
        start = content.index(START_EXTERNAL) + 1
        end = content.index(END_EXTERNAL)
        return content[start:end]

    imports = []
    for m in modules:
        imports.extend(external_modules(m))
    # Remove duplicates and sort
    imports = list(set(imports))
    imports.sort()
    return imports

def remove_module_imports(content):
    start = content.index(START_EXTERNAL)
    end = content.index(END_EXTERNAL) + 1
    for _ in range(start, end):
        content.pop(start)
    start = content.index(START_INTERNAL)
    end = content.index(END_INTERNAL) + 1
    for _ in range(start, end):
        content.pop(start)

def remove_imports(modules):
    for m in modules:
        remove_module_imports(m)

def parse_code_or_comment(content):
    comment = ''
    # Parse comment above code
    while content:
        line = content[0]
        if not (re.match('^( )*#', line) or line.isspace()):
            break
        comment += line
        content.pop(0)
    code = ''
    first_run = True
    # Parse code
    while content:
        line = content[0]
        if line.startswith('"""') or line.startswith('#'):
            break
        if not first_run and (line.startswith('def ') or line.startswith('class')):
            break
        code += line
        first_run = False
        content.pop(0)
    nb['cells'].append(nbf.v4.new_code_cell((comment + code).rstrip()))

def parse_markdown(content):
    markdown = content[0][3:]  # Text immediately after """
    content.pop(0)  # remove beginning """
    while content:
        line = content.pop(0)
        if line.startswith('"""'):
            nb['cells'].append(nbf.v4.new_markdown_cell(markdown))
            break
        markdown += line

def parse(content):
    while content:
        line = content[0]
        if line.isspace():
            # Skip whitespace.
            content.pop(0)
        elif line.startswith('# python_modules_to_import'):
            # interrupt to include modules
            content.pop(0)
            break
        elif line.startswith('"""'):
            parse_markdown(content)
        else:
            parse_code_or_comment(content)

with open('../nanodip.py', 'r') as f:
    nanodip = f.readlines()
with open('../config.py', 'r') as f:
    config = f.readlines()
with open('../utils.py', 'r') as f:
    utils = f.readlines()
with open('../data.py', 'r') as f:
    data = f.readlines()
with open('../plots.py', 'r') as f:
    plots = f.readlines()

nanodip_modules = [nanodip, config, utils, data, plots]
imports = external_imports(nanodip_modules)
remove_imports(nanodip_modules)

nb = nbf.v4.new_notebook()

parse(nanodip)
for m in [imports] + nanodip_modules[1:] + [nanodip]:
    parse(m)

nbf.write(nb, 'nanodip.ipynb')
copy_tree('../static','static')
copy_tree('../templates','templates')
