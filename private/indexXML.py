#!/usr/bin/env python
# Copyright (C) 2017 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.
# Author: Philip Ilten, May 2017.

# This script is used to build the search index, excute as ./indexXML.py.
import glob, re
import xml.etree.cElementTree as xml

# Clean the special characters in a line and append to a list.
def append(text, line):
    if not line: return
    for f, r in [('\n', ''), ('@', '\@'), ('AMPERSAND', '&'), ('&lt;', '<'),
                 ('&gt;', '>'), ('"', '\'')]:
        line = line.replace(f, r)
    line = line.strip()
    if not line: return
    text += [line]

# Parse an XML element and write it to the index.
def parse(index, base, element, key = None, lnk = 0):
    if element.tag in ['particle', 'channel']: return lnk
    if 'name' in element.attrib and element.tag not in ['argument']:
        key = element.attrib['name']
        if element.tag not in ['chapter']: lnk += 1
    if key not in index:
        index[key] = {'link':'%s.html' % base, 'text': []}
        if element.tag != 'chapter': index[key]['link'] += '#anchor%i' % lnk
    append(index[key]['text'], element.text)
    for child in element: lnk = parse(index, base, child, key, lnk)
    append(index[key]['text'], element.tail)
    return lnk

# The main script.
names  = glob.glob('../share/Pythia8/xmldoc/*.xml')
names  = sorted(names)
ignore = ['Bibliography', 'Frontpage', 'Glossary', 'Index', 'ProgramClasses',
          'ProgramFiles', 'ProgramMethods', 'SaveSettings', 'UpdateHistory',
          'Version', 'Welcome']
index  = {}
for name in names:
    xmldoc = file(name).read().replace('&', 'AMPERSAND')
    name = name.split('/')[-1].replace('.xml', '')
    if name in ignore: continue
    tree = xml.fromstring(xmldoc)
    parse(index, name, tree)
out = file('../share/Pythia8/htmldoc/Index.js', 'w')
out.write('var index = [')
for idx, (key, val) in enumerate(index.iteritems()):
    out.write('%s{"name":"%s","link":"%s","text":"%s"}' % (
            (',' if idx else ''), key, val['link'], 
            re.sub(r'\s([ ?.!"](?:\s|$))', r'\1', ' '.join(val['text']))))
out.write('];')
out.close()
