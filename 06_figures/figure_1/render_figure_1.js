// Renders Figure_1.svg to a PNG of 2067 x 2766 px (175 mm wide at 300 dpi) and to a vector PDF of 175 mm x 234.167 mm.
// Usage: node render_figure_1.js Figure_1.svg OUT_PREFIX
// Needs Node.js with the package "playwright" and its Chromium. Writes OUT_PREFIX.png and OUT_PREFIX.pdf; existing files are not overwritten.
const { chromium } = require('playwright');
const fs = require('fs');
(async () => {
  const [svgPath, out] = process.argv.slice(2);
  if (!svgPath || !out) { console.error('usage: node render_figure_1.js Figure_1.svg OUT_PREFIX'); process.exit(1); }
  for (const ext of ['.png', '.pdf']) {
    if (fs.existsSync(out + ext)) { console.error('STOP: output exists: ' + out + ext); process.exit(1); }
  }
  const svg = fs.readFileSync(svgPath, 'utf8');
  const src = 'data:image/svg+xml;base64,' + Buffer.from(svg, 'utf8').toString('base64');
  const browser = await chromium.launch();
  const W = 2067, H = 2766;
  const page = await browser.newPage({ viewport: { width: W, height: H }, deviceScaleFactor: 1 });
  await page.setContent(`<html><body style="margin:0;background:white"><img id=i style="width:${W}px;height:${H}px;display:block"></body></html>`);
  await page.evaluate(s => { document.getElementById('i').src = s; }, src);
  await page.waitForTimeout(1500);
  await page.screenshot({ path: out + '.png' });
  const sheet = await browser.newPage();
  await sheet.setContent(`<html><head><style>@page{size:175mm 234.167mm;margin:0} html,body{margin:0;padding:0;background:white} img{display:block;width:175mm;height:234.167mm}</style></head><body><img id=i></body></html>`);
  await sheet.evaluate(s => { document.getElementById('i').src = s; }, src);
  await sheet.waitForTimeout(1500);
  await sheet.pdf({ path: out + '.pdf', width: '175mm', height: '234.167mm', printBackground: true, pageRanges: '1', margin: { top: '0', right: '0', bottom: '0', left: '0' } });
  await browser.close();
})();
