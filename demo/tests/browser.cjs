/* Real exported results only. No simulated scientific output is supplied by this test. */
const { chromium } = require('playwright');
const assert = require('node:assert/strict');
const path = require('node:path');
const fs = require('node:fs/promises');

(async () => {
  const output = process.env.DEMO_SCREENSHOTS || path.resolve('test-results');
  await fs.mkdir(output, { recursive: true });
  const browser = await chromium.launch({ headless: true, ...(process.env.PLAYWRIGHT_CHANNEL ? {channel:process.env.PLAYWRIGHT_CHANNEL} : {}), args: ['--enable-webgl', '--use-gl=angle', '--use-angle=swiftshader', '--no-sandbox'] });
  const errors = [];
  const checks = [];
  for (const [name, width, height] of [['desktop',1440,1000],['mobile',390,844]]) {
    const page = await browser.newPage({ viewport:{width,height}, deviceScaleFactor:1 });
    page.on('pageerror',e=>errors.push(e.message));
    await page.goto(process.env.DEMO_URL || 'http://127.0.0.1:5173/');
    await page.waitForSelector('canvas[data-residues="444"]');
    await page.waitForTimeout(700);
    const pixels = await page.locator('canvas').evaluate(canvas => {
      const copy = document.createElement('canvas'); copy.width=canvas.width;copy.height=canvas.height;
      const context=copy.getContext('2d');context.drawImage(canvas,0,0);
      const pixels=context.getImageData(0,0,copy.width,copy.height).data;
      let protein=0,ligand=0;
      for(let i=0;i<pixels.length;i+=4){if(pixels[i+1]>pixels[i]*1.3&&pixels[i+1]>pixels[i+2]*1.04)protein++;if(pixels[i]>pixels[i+1]*1.3&&pixels[i]>pixels[i+2]*1.25)ligand++;}
      return {protein,ligand,total:copy.width*copy.height};
    });
    assert(pixels.protein>500, `${name}: protein pixels absent`);
    assert(pixels.ligand>10, `${name}: ligand pixels absent`);
    const before = await page.locator('canvas').screenshot();
    await page.getByRole('button',{name:'Rotate structure',exact:true}).click();
    await page.waitForTimeout(700);
    const after=await page.locator('canvas').screenshot();
    assert(!before.equals(after),`${name}: rotation did not change canvas`);
    await page.getByRole('button',{name:'Rotate structure',exact:true}).click();
    await page.getByRole('button',{name:'Overlay another aligned seed'}).click();
    await page.waitForTimeout(200);
    assert.equal(await page.locator('#overlay-key').isVisible(),true);
    await page.getByRole('button',{name:'Overlay another aligned seed'}).click();
    await page.getByRole('button',{name:'Focus ligand pocket'}).click();
    const pocket=await page.locator('canvas').screenshot();
    assert(!pocket.equals(before),`${name}: pocket focus had no effect`);
    await page.getByRole('button',{name:'Fit full structure'}).click();
    await page.screenshot({path:path.join(output,`${name}-structure.png`),fullPage:true});
    assert(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth),`${name}: horizontal page overflow`);
    await page.locator('#subject').selectOption('ecoli_G96A');
    assert.equal(await page.locator('canvas').getAttribute('data-residues'),'427');
    await page.locator('#context').selectOption('s3p_only');
    await page.locator('#seed').selectOption('107');
    assert.match(await page.locator('canvas').getAttribute('aria-label'),/s3p_only, seed 107/);
    await page.getByRole('button',{name:'Binding calibration',exact:true}).click();
    await page.getByText('Direction warning',{exact:true}).waitFor();
    await page.locator('.scientific-figure').waitFor();
    assert.equal(await page.locator('.scientific-figure').evaluate(img=>img.complete&&img.naturalWidth>0),true);
    await page.screenshot({path:path.join(output,`${name}-calibration.png`),fullPage:true});
    await page.getByRole('button',{name:'Candidate archive',exact:true}).click();
    assert.equal(await page.locator('tbody tr').count(),7);
    await page.getByRole('button',{name:'Evidence & review',exact:true}).click();
    await page.getByText('API request verified',{exact:true}).waitFor();
    await page.screenshot({path:path.join(output,`${name}-evidence.png`),fullPage:true});
    await page.locator('[data-stage="2"]').click();
    assert.equal(await page.locator('.stage-detail').isVisible(),true);
    await page.getByRole('button',{name:'Close stage details'}).click();
    const other=await page.locator('#target option').evaluateAll(options=>options.find(o=>o.value!=='AT2G45300').value);
    await page.locator('#target').selectOption(other);
    await page.getByRole('heading',{name:'No run available'}).waitFor();
    await page.getByRole('button',{name:'Return to EPSPS'}).click();
    await page.getByRole('button',{name:'Structure & results',exact:true}).click();
    await page.waitForSelector('canvas');
    checks.push({viewport:name,pixels,interactions:'passed'});
    await page.close();
  }
  assert.deepEqual(errors,[]);
  await browser.close();
  await fs.writeFile(path.join(output,'browser-checks.json'),JSON.stringify({checks,errors},null,2));
  console.log(JSON.stringify({checks,errors},null,2));
})().catch(error=>{console.error(error);process.exit(1);});
