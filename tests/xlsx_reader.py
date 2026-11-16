"""
MODULE:     t.graminr (tests)

AUTHOR(S):  Yann Chemin

PURPOSE:    Minimal reader for the GrasSim reference workbooks shipped with
            GraMinR. openpyxl cannot open them (they carry a drawing
            relationship pointing at a part that is not in the archive), and
            only the cached cell values are needed, so the sheet XML is read
            directly.

SPDX-License-Identifier: Unlicense
"""

import xml.etree.ElementTree as ET
import zipfile

MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
NS = {"m": MAIN}


def _column(reference):
    """Column letters of a cell reference such as "BC12"."""
    return "".join(char for char in reference if char.isalpha())


def _column_index(letters):
    value = 0
    for char in letters:
        value = value * 26 + (ord(char) - ord("A") + 1)
    return value


def read_sheet(path, sheet_name):
    """Rows of one worksheet as a list of dicts keyed by the header row."""
    with zipfile.ZipFile(path) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relations = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {}
        for relation in relations:
            if relation.get("Type", "").endswith("/worksheet"):
                targets[relation.get("Id")] = relation.get("Target")
        part = None
        for sheet in workbook.iter("{%s}sheet" % MAIN):
            rid = sheet.get(
                "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
            )
            if sheet.get("name") == sheet_name:
                part = targets[rid]
        if part is None:
            raise KeyError("no sheet named {!r} in {}".format(sheet_name, path))

        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            table = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in table.findall("m:si", NS):
                shared.append(
                    "".join(node.text or "" for node in item.iter("{%s}t" % MAIN))
                )

        sheet_xml = ET.fromstring(archive.read("xl/" + part.lstrip("/")))

    rows = []
    for row in sheet_xml.find("m:sheetData", NS):
        cells = {}
        for cell in row.findall("m:c", NS):
            value = cell.find("m:v", NS)
            if value is None:
                inline = cell.find("m:is", NS)
                text = (
                    "".join(node.text or "" for node in inline.iter("{%s}t" % MAIN))
                    if inline is not None
                    else None
                )
            elif cell.get("t") == "s":
                text = shared[int(value.text)]
            else:
                text = value.text
            cells[_column(cell.get("r"))] = text
        rows.append(cells)

    if not rows:
        return []
    header = rows[0]
    order = sorted(header, key=_column_index)
    names = [header[letters] for letters in order]
    return [
        {names[index]: row.get(order[index]) for index in range(len(order))}
        for row in rows[1:]
    ]
