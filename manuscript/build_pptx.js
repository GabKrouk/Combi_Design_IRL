// Figures deck: one A4 portrait page per figure (title, figure at full width, legend underneath).
const pptxgen = require('pptxgenjs'); const fs = require('fs'); const path = require('path');
const {execSync} = require('child_process');
const FIG = process.env.FIG_DIR || '../analysis/paper';
const D = JSON.parse(fs.readFileSync('ms.json', 'utf8'));

const p = new pptxgen();
p.defineLayout({name: 'A4_PORTRAIT', width: 8.27, height: 11.69});
p.layout = 'A4_PORTRAIT';
p.title = 'Figures - Krouk & Shasha'; p.author = 'Gabriel Krouk, Dennis E. Shasha';
const PW = 8.27, PH = 11.69, MX = 0.5, FONT = 'Arial';

// rich text with *italic* / **bold** / ^sup^ / ~sub~
function rich(t, base = {}) {
  const out = []; const re = /(\*\*[^*]+?\*\*|\*[^*]+?\*|\^[^^]+?\^|~[^~]+?~)/g; let last = 0, m;
  const push = (text, o) => text && out.push({text, options: {...base, ...o}});
  while ((m = re.exec(t)) !== null) {
    push(t.slice(last, m.index), {});
    const s = m[0];
    if (s.startsWith('**')) rich(s.slice(2, -2), {...base, bold: true}).forEach(x => out.push(x));
    else if (s.startsWith('*')) rich(s.slice(1, -1), {...base, italic: true}).forEach(x => out.push(x));
    else if (s.startsWith('^')) push(s.slice(1, -1).replace(/\*/g, ''), {superscript: true});
    else push(s.slice(1, -1).replace(/\*/g, ''), {subscript: true});
    last = m.index + s.length;
  }
  push(t.slice(last), {});
  return out;
}
const plain = t => t.replace(/\*\*|\*|\^|~/g, '');

// cover page
let s = p.addSlide(); s.background = {color: 'FFFFFF'};
s.addText(plain(D.meta.TITLE), {x: MX, y: 3.2, w: PW - 2 * MX, h: 1.6, fontFace: FONT, fontSize: 24, bold: true, color: '1F6F5C', valign: 'top'});
s.addText('Gabriel Krouk and Dennis E. Shasha', {x: MX, y: 5.0, w: PW - 2 * MX, h: 0.4, fontFace: FONT, fontSize: 14, color: '333333'});
s.addText('Figures 1–5 (main text) and Figures S1–S3 (SI Appendix)', {x: MX, y: 5.5, w: PW - 2 * MX, h: 0.4, fontFace: FONT, fontSize: 12, color: '666666'});

const L = D.legends;
for (let i = 0; i < L.length; i++) {
  if (L[i].t !== 'h2') continue;
  const m = L[i].x.match(/^Fig\. (S?\d+)\./); const f = path.join(FIG, `Fig${m[1]}.png`);
  const body = (L[i + 1] && L[i + 1].t === 'p') ? L[i + 1].x : '';
  const [w, h] = execSync(`python3 -c "from PIL import Image; print(*Image.open('${f}').size)"`).toString().trim().split(' ').map(Number);
  s = p.addSlide(); s.background = {color: 'FFFFFF'};
  s.addText(rich(L[i].x, {fontFace: FONT, fontSize: 13, bold: true, color: '222222'}), {x: MX, y: 0.4, w: PW - 2 * MX, h: 0.5, valign: 'top', margin: 0});
  const W = PW - 2 * MX, H = W * h / w, y0 = 1.0;
  s.addImage({path: f, x: MX, y: y0, w: W, h: H});
  const yL = y0 + H + 0.25;
  s.addText(rich(body, {fontFace: FONT, fontSize: 9, color: '333333'}), {x: MX, y: yL, w: W, h: PH - yL - 0.5, valign: 'top', margin: 0, paraSpaceAfter: 0, lineSpacingMultiple: 1.1});
  s.addText(String(m[1]), {x: PW - MX - 1, y: PH - 0.45, w: 1, h: 0.3, fontFace: FONT, fontSize: 8, color: '999999', align: 'right'});
  s.addNotes(plain(L[i].x) + '\n' + plain(body));
}
p.writeFile({fileName: 'Krouk_Shasha_figures_A4.pptx'}).then(() => console.log('ok'));
