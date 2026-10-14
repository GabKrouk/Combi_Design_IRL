const fs = require('fs');
const {Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType, Footer, PageNumber, ImageRun,
  Table, TableRow, TableCell, WidthType, BorderStyle, ShadingType, LineNumberRestartFormat, PageBreak} = require('docx');

const D = JSON.parse(fs.readFileSync('ms.json', 'utf8'));
const FONT = 'Arial';
const SP = {line: 480};                         // double spacing for review

// inline markup: **bold**, *italic*, ^sup^, ~sub~ (recursive)
function runs(text, st = {}) {
  const out = []; const re = /(\*\*[^*]+?\*\*|\*[^*]+?\*|\^[^^]+?\^|~[^~]+?~)/g; let last = 0, m;
  while ((m = re.exec(text)) !== null) {
    if (m.index > last) out.push(new TextRun({text: text.slice(last, m.index), font: FONT, ...st}));
    const t = m[0];
    if (t.startsWith('**')) out.push(...runs(t.slice(2, -2), {...st, bold: true}));
    else if (t.startsWith('*')) out.push(...runs(t.slice(1, -1), {...st, italics: true}));
    else if (t.startsWith('^')) out.push(...runs(t.slice(1, -1), {...st, superScript: true}));
    else out.push(...runs(t.slice(1, -1), {...st, subScript: true}));
    last = m.index + t.length;
  }
  if (last < text.length) out.push(new TextRun({text: text.slice(last), font: FONT, ...st}));
  return out;
}
const P = (text, o = {}) => new Paragraph({children: runs(text, o.run || {}), spacing: {...SP, after: o.after ?? 0},
  indent: o.indent, alignment: o.align, keepNext: o.keepNext});
const H1 = t => new Paragraph({heading: HeadingLevel.HEADING_1, children: [new TextRun({text: t, font: FONT, bold: true, size: 26})],
  spacing: {before: 360, after: 120, line: 360}, keepNext: true});
const H2 = t => new Paragraph({heading: HeadingLevel.HEADING_2, children: runs(t, {bold: true, size: 23}),
  spacing: {before: 240, after: 60, line: 360}, keepNext: true});

const kids = [];
const M = D.meta;
kids.push(new Paragraph({children: runs(M.TITLE, {bold: true, size: 32}), spacing: {after: 240, line: 360}}));
kids.push(P(M.AUTHORS, {run: {size: 24}, after: 120}));
kids.push(P(M.AFFIL, {run: {size: 20}, after: 120}));
kids.push(P(M.CORR, {run: {size: 20}, after: 120}));
kids.push(P('**Classification:** ' + M.CLASS, {run: {size: 20}}));
kids.push(P('**Keywords:** ' + M.KEYWORDS, {run: {size: 20}}));
kids.push(P('**This PDF/Word file includes:** Main text, Figures 1 to 5, Table 1. **Supporting Information:** Figures S1 to S3 (legends at the end of this file).', {run: {size: 20}}));

let first = true;
for (const b of D.blocks) {
  if (b.t === 'h1') {
    if (['Abstract', 'Significance', 'Introduction'].includes(b.x)) kids.push(new Paragraph({children: [new PageBreak()]}));
    kids.push(H1(b.x));
  } else if (b.t === 'h2') kids.push(H2(b.x));
  else {
    const isEq = b.x.startsWith('*N*~95~ ≈');
    kids.push(P(b.x, {after: 120, align: isEq ? AlignmentType.CENTER : undefined}));
  }
}

// references
kids.push(new Paragraph({children: [new PageBreak()]}));
kids.push(H1('References'));
D.refs.forEach((r, i) => kids.push(new Paragraph({children: runs(`${i + 1}. ` + r, {size: 20}), spacing: {line: 360, after: 60},
  indent: {left: 440, hanging: 440}})));

// table 1
kids.push(new Paragraph({children: [new PageBreak()]}));
kids.push(new Paragraph({children: runs('**' + D.table.title + '**', {size: 20}), spacing: {after: 120, line: 300}}));
const W = [1500, 1700, 1500, 1500, 1600, 1560]; const TW = W.reduce((a, b) => a + b, 0);
const border = {style: BorderStyle.SINGLE, size: 4, color: '999999'};
const none = {style: BorderStyle.NONE, size: 0, color: 'FFFFFF'};
const trows = D.table.rows.map((r, i) => new TableRow({tableHeader: i === 0, children: r.map((c, j) => new TableCell({
  width: {size: W[j], type: WidthType.DXA}, margins: {top: 60, bottom: 60, left: 80, right: 80},
  shading: i === 0 ? {fill: 'E8F1EE', type: ShadingType.CLEAR, color: 'auto'} : undefined,
  borders: {top: i === 0 ? border : none, bottom: (i === 0 || i === D.table.rows.length - 1) ? border : none, left: none, right: none},
  children: [new Paragraph({alignment: j === 0 ? AlignmentType.LEFT : AlignmentType.CENTER,
    children: [new TextRun({text: c, font: FONT, size: 19, bold: i === 0 || j === 5, color: j === 5 && i > 0 ? '1F6F5C' : undefined})]})]}))}));
kids.push(new Table({width: {size: TW, type: WidthType.DXA}, columnWidths: W, rows: trows}));
kids.push(new Paragraph({children: runs(D.table.note, {size: 18}), spacing: {before: 120, line: 300}}));

// figures with legends
const sizes = {Fig1: [6.5, 2.63], Fig2: [6.5, 4.31], Fig3: [6.5, 4.93], Fig4: [6.5, 4.2], Fig5: [6.5, 2.26]};
const img = require('child_process');
for (const b of D.legends) {
  if (b.t === 'h2') {
    const m = b.x.match(/^Fig\. (S?\d+)\./); const name = 'Fig' + m[1];
    kids.push(new Paragraph({children: [new PageBreak()]}));
    if (m[1].startsWith('S') && name === 'FigS1') kids.push(H1('Supporting Information figure legends (SI Appendix)'));
    const f = `${process.env.FIG_DIR || '../analysis/paper'}/${name}.png`;
    if (fs.existsSync(f)) {
      const dim = img.execSync(`python3 -c "from PIL import Image; im=Image.open('${f}'); print(im.size[0], im.size[1])"`).toString().trim().split(' ').map(Number);
      const wpx = 620, hpx = Math.round(wpx * dim[1] / dim[0]);
      kids.push(new Paragraph({alignment: AlignmentType.CENTER, children: [new ImageRun({type: 'png', data: fs.readFileSync(f),
        transformation: {width: wpx, height: hpx}, altText: {title: name, description: b.x, name}})]}));
    }
    kids.push(new Paragraph({children: runs('**' + b.x + '**', {size: 20}), spacing: {before: 160, after: 60, line: 300}, keepNext: true}));
  } else kids.push(new Paragraph({children: runs(b.x, {size: 20}), spacing: {line: 300, after: 60}}));
}

const doc = new Document({
  creator: 'Gabriel Krouk; Dennis E. Shasha', title: M.TITLE,
  styles: {default: {document: {run: {font: FONT, size: 22}}},
    paragraphStyles: [
      {id: 'Heading1', name: 'Heading 1', basedOn: 'Normal', next: 'Normal', quickFormat: true, run: {font: FONT, size: 26, bold: true}, paragraph: {outlineLevel: 0}},
      {id: 'Heading2', name: 'Heading 2', basedOn: 'Normal', next: 'Normal', quickFormat: true, run: {font: FONT, size: 23, bold: true}, paragraph: {outlineLevel: 1}}]},
  sections: [{
    properties: {page: {size: {width: 11906, height: 16838}, margin: {top: 1418, bottom: 1418, left: 1418, right: 1418}},
      lineNumbers: {countBy: 1, restart: LineNumberRestartFormat.CONTINUOUS, distance: 360}},
    footers: {default: new Footer({children: [new Paragraph({alignment: AlignmentType.CENTER,
      children: [new TextRun({children: [PageNumber.CURRENT], font: FONT, size: 18})]})]})},
    children: kids}]});
Packer.toBuffer(doc).then(b => { fs.writeFileSync('Krouk_Shasha_manuscript.docx', b); console.log('ok', b.length); });
