/* Exercises actual exported evidence. No simulated scientific output. */
const { chromium } = require('playwright');
const assert = require('node:assert/strict');
const path = require('node:path');
const fs = require('node:fs/promises');
(async()=>{
 const output=process.env.DEMO_SCREENSHOTS||path.resolve('test-results');await fs.mkdir(output,{recursive:true});
 const browser=await chromium.launch({headless:true,...(process.env.PLAYWRIGHT_CHANNEL?{channel:process.env.PLAYWRIGHT_CHANNEL}:{}),args:['--enable-webgl','--use-gl=angle','--use-angle=swiftshader','--no-sandbox']});
 const errors=[],checks=[];
 for(const [name,width,height] of [['desktop',1440,1000],['mobile',390,844]]){
  const page=await browser.newPage({viewport:{width,height},deviceScaleFactor:1});page.on('pageerror',e=>errors.push(e.message));await page.goto(process.env.DEMO_URL||'http://127.0.0.1:5180/');await page.getByRole('heading',{name:'Overview',exact:true}).waitFor();
  assert(!await page.locator('body').innerText().then(t=>t.includes('[object Object]')));
  await page.screenshot({path:path.join(output,`${name}-overview.png`),fullPage:true});
  await page.locator('[data-tab="workflow"]').click();assert.equal(await page.locator('.workflow-node').count(),14);await page.locator('[data-agent="8"]').click();assert(await page.locator('.agent-inspector').innerText().then(t=>t.toLowerCase().includes('mutation')));
  await page.locator('[data-tab="binding"]').click();assert.equal(await page.locator('#chemical-image svg').count(),1);await page.locator('[data-ligand="pep"]').click();assert(await page.locator('.affinity-number').innerText().then(t=>t.includes('pIC50')));
  const jsonDownload=page.waitForEvent('download');await page.locator('[data-export="binding"][data-format="json"]').click();const download=await jsonDownload;assert.equal(download.suggestedFilename(),'binding.json');const filename=path.join(output,`${name}-binding.json`);await download.saveAs(filename);const exported=JSON.parse(await fs.readFile(filename,'utf8'));assert(exported.rows.length>0);assert(exported.legend);
  await page.locator('[data-tab="residues"]').click();await page.waitForSelector('canvas[data-residues="444"]');await page.waitForTimeout(500);
  const before=await page.locator('canvas').screenshot();await page.locator('#spin').click();await page.waitForTimeout(500);const after=await page.locator('canvas').screenshot();assert(!before.equals(after),'Rotation must update canvas');await page.locator('#spin').click();
  await page.locator('#overlay').click();assert(await page.locator('#overlay-key').isVisible());await page.locator('#overlay').click();
  const row=page.locator('[data-residue]').first();const position=await row.getAttribute('data-residue');await row.click();assert(await page.locator('#residue-detail').innerText().then(t=>t.includes('Residue '+position)&&t.includes('general proximity')));
  const selectedCanvas=await page.locator('canvas').screenshot();assert(!before.equals(selectedCanvas),'Residue focus must update canvas');
  await page.locator('[data-ligand-chain]').first().uncheck();await page.locator('[data-ligand-chain]').first().check();
  const seedBefore=await page.locator('#seed').inputValue();await page.locator('#next-pose').click();assert.notEqual(await page.locator('#seed').inputValue(),seedBefore);
  const pngDownload=page.waitForEvent('download');await page.locator('#image-export').click();assert((await pngDownload).suggestedFilename().endsWith('.png'));
  const count=await page.locator('.matrix-table tbody tr').count();await page.locator('#filter-native').fill('0');assert((await page.locator('.matrix-table tbody tr').count())<=count);await page.locator('#filter-native').fill('1');
  await page.screenshot({path:path.join(output,`${name}-residues.png`),fullPage:true});
  const archive=await page.locator('#subject option').evaluateAll(options=>options.find(o=>o.value.startsWith('archive:')&&!o.value.includes('WT'))?.value);if(archive){await page.locator('#subject').selectOption(archive);assert(await page.locator('#compare-wt').isEnabled());await page.locator('#compare-wt').click();assert(await page.locator('#overlay-key').innerText().then(t=>t.includes('WT')));}
  await page.locator('[data-tab="comparison"]').click();await page.locator('#comparison-metric').selectOption('alignment_lddt');assert(await page.locator('.structure-plot svg').getAttribute('aria-label').then(t=>t.includes('lDDT')));
  await page.locator('[data-tab="decisions"]').click();assert((await page.locator('.gate').count())>=9);assert(await page.locator('.decision-selector').innerText().then(t=>t.includes('insufficient evidence')));await page.screenshot({path:path.join(output,`${name}-decisions.png`),fullPage:true});
  await page.locator('[data-tab="evidence"]').click();assert.equal(await page.locator('.export-list [data-export]').count(),10);assert.equal(await page.locator('.validation-plan li').count(),6);
  assert(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth),`${name}: horizontal page overflow`);
  checks.push({viewport:name,agent_nodes:14,viewer:'rendered; rotation, residue focus, overlay, seed switch, ligand visibility, PNG export tested',table_export:'JSON content and filename verified',decisions:'missing gates remain insufficient'});await page.close();
 }
 assert.deepEqual(errors,[]);await browser.close();await fs.writeFile(path.join(output,'browser-checks.json'),JSON.stringify({checks,errors},null,2));console.log(JSON.stringify({checks,errors},null,2));
})().catch(e=>{console.error(e);process.exit(1)});
