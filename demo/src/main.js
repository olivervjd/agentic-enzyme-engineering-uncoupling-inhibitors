import { createIcons, Sprout, LockKeyhole, Download, X, Scan, Focus, Rotate3d, Layers2,
  Box, ChartNoAxesCombined, Archive, FileSearch, History, ArrowUpRight, CircleAlert,
  GitBranch, BookOpen, ScanLine, Layers, Network, Dna, Scale, Files, NotebookPen, BadgeCheck } from 'lucide';
import { createViewer } from './viewer.js';
import './style.css';

const root = document.querySelector('#app');
const icons = { Sprout, LockKeyhole, Download, X, Scan, Focus, Rotate3d, Layers2,
  Box, ChartNoAxesCombined, Archive, FileSearch, History, ArrowUpRight, CircleAlert,
  GitBranch, BookOpen, ScanLine, Layers, Network, Dna, Scale, Files, NotebookPen, BadgeCheck };
const esc = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const icon = name => `<i data-lucide="${name}" aria-hidden="true"></i>`;
const fmt = (x, places=3) => Number.isFinite(x) ? x.toFixed(places) : 'Not measured';
const pretty = x => String(x).replaceAll('_', ' ').toLowerCase();
let data, structures, viewer;
let tab = 'structure', subject = 'arabidopsis_wt', context = 'herbicide', seed = 101, stage = null;
let target = 'AT2G45300', spinning = false, overlay = false;
const agentDefs = [
  ['Workflow Orchestrator','workflow_orchestrator','git-branch','Coordinates evidence and holds the run when required checks are unresolved.'],
  ['Evidence Agent','evidence','book-open','Literature retrieval and model synthesis have separate statuses. Retrieved abstracts are not experimental validation.'],
  ['Structure-Quality Agent','structure_quality','scan-line','Structure comparisons pass; chemistry, target-specific accuracy and functional evidence still need review.'],
  ['Pose-Ensemble Agent','pose_ensemble','layers','Boltz cofolded poses are inspected independently of candidate approval. Independent DiffDock availability is reported separately.'],
  ['Interaction-Fingerprint Agent','interaction_fingerprint','network','Compares herbicide and native-ligand contacts. Shared and native-function sites remain protected.'],
  ['Constrained Mutation Agent','constrained_mutation','dna','No new mutations proposed in this run. The archive contains six earlier candidates.'],
  ['Multi-Oracle Scoring Agent','candidate_oracles','scale','Requires mutation-specific structures, matched ligand contexts and valid affinity evidence.'],
  ['Review Packet Agent','workflow_review_packet','files','Diagnostic packets capture failures and missing evidence even when candidate design is blocked. Candidate reviews remain separate.'],
  ['Learning Agent','learning','notebook-pen','Awaiting experimental assay results. No learning update has occurred.'],
  ['LLM Judge / Evaluation','workflow_llm_judge','badge-check','The workflow judge audits the diagnostic packet and a separate reviewer response. This is not candidate approval.']
];
function stageRecord(name) {
  const aliases={pose_ensemble:'pose_ensemble_and_fingerprint',interaction_fingerprint:'pose_ensemble_and_fingerprint',workflow_review_packet:'scoring_and_review_packets',workflow_llm_judge:'llm_judge'};
  return data.workflow.find(r=>r.stage===name) || data.workflow.find(r=>r.stage===aliases[name]);
}
function status(name) { return name==='workflow_orchestrator' ? data.run.status : stageRecord(name)?.status || 'NOT_RUN'; }
function badge(text, tone='neutral') { return `<span class="badge ${tone}">${esc(text)}</span>`; }
function iconButton(id, name, label, pressed=null) {
  return `<button id="${id}" class="icon-button" title="${label}" aria-label="${label}"${pressed === null ? '' : ` aria-pressed="${pressed}"`}>${icon(name)}</button>`;
}
function legend(title, text, open=false) { return `<details class="legend"${open ? ' open' : ''}><summary>${esc(title)}</summary><p>${esc(text)}</p></details>`; }
function download(name, content) {
  const url = URL.createObjectURL(new Blob([content], {type:'application/json'}));
  const a = document.createElement('a'); a.href = url; a.download = name; a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
function shell() {
  const available = target === 'AT2G45300';
  root.innerHTML = `<a class="skip-link" href="#main">Skip to results</a>
    <header class="app-header"><div class="brand">${icon('sprout')}<span>Resistance <b>Workflow</b></span></div>
      <div class="header-right"><span class="recorded"><span></span>Recorded scientific run</span>${iconButton('download','download','Download displayed results')}</div></header>
    <div class="layout"><aside class="sidebar"><label for="target">TARGET</label><select id="target">${data.targets.map(t => `<option value="${esc(t.id)}" ${target===t.id?'selected':''}>${esc(t.id)} / ${esc(t.herbicide)}</option>`).join('')}</select>
      <div class="rail-heading"><span>AGENT WORKFLOW</span><span>10 stages</span></div>
      <nav class="agents" aria-label="Workflow agents">${agentDefs.map(([name,key,symbol], i) => {
        const value = available ? status(key) : 'NOT_CONFIGURED';
        const tone = value === 'COMPLETED' ? 'complete' : value === 'NEEDS_REASSESSMENT' ? 'attention' : 'pending';
        return `<button class="agent ${tone}${stage===i?' selected':''}" data-stage="${i}" ${available?'':'disabled'}><span class="stage-number">${String(i+1).padStart(2,'0')}</span><span class="agent-text">${esc(name)}<small>${esc(pretty(value))}</small></span><span class="stage-dot"></span></button>`;
      }).join('')}</nav><div class="rail-footer">${icon('lock-keyhole')}<span>Review gates active<br><strong>No automatic candidate approval</strong></span></div></aside>
      <main id="main"><div class="page-title"><div><div class="eyebrow">ARABIDOPSIS / COMPUTATIONAL DESIGN</div><h1>${available?'EPSPS <span>/</span> Glyphosate':esc(data.targets.find(t=>t.id===target)?.name || target)}</h1><p>${available?'AT2G45300 · Mature protein 77–520 · Matched-MSA calibration':'No recorded scientific run for this target'}</p></div>${badge(available?'Needs reassessment':'Not configured','warning')}</div>
      ${available ? `<div class="metrics"><div><span>Structure predictions</span><strong>${data.run.predictions}<small>Boltz-2 · three seeds</small></strong></div><div><span>Structure checks passed</span><strong>${data.run.passed}<em>/ ${data.run.comparisons}</em><small>Repeatability + homolog controls</small></strong></div><div><span>Experimental binding calibration</span><strong class="metric-status">${data.binding.status==='VALIDATED_WITHIN_SCOPE'?'Scope-limited':'Not validated'}<small>${data.binding.eligible} eligible matched observations</small></strong></div><div><span>Candidates advanced</span><strong>0<small>Scientific gates remain closed</small></strong></div></div>
      <nav class="tabs" aria-label="Results views">${[['structure','Structure & results','box'],['calibration','Binding calibration','chart-no-axes-combined'],['archive','Candidate archive','archive'],['evidence','Evidence & review','file-search']].map(([id,label,symbol])=>`<button data-tab="${id}" aria-current="${tab===id?'page':'false'}">${icon(symbol)}${label}</button>`).join('')}</nav>
      ${stage !== null ? `<section class="stage-detail"><div>${icon(agentDefs[stage][2])}<strong>${esc(agentDefs[stage][0])}</strong>${badge(pretty(status(agentDefs[stage][1])))}</div><p>${esc(agentDefs[stage][3])}</p>${(stageRecord(agentDefs[stage][1])?.reasons||[]).map(r=>`<p>${esc(r)}</p>`).join('')}${iconButton('close-stage','x','Close stage details')}</section>` : ''}
      <div id="view"></div>` : `<section class="empty"><h2>No run available</h2><p>Only EPSPS has real recorded predictions in this demo. No results have been generated for ${esc(target)}.</p><button id="return-epsps" class="primary">Return to EPSPS</button></section>`}
      <footer class="page-footer"><span>Recorded run: ${esc(data.run.id)}</span><span>Predictions are not experimental proof.</span></footer></main></div>`;
  document.querySelector('#download').onclick = () => download('epsps-displayed-results.json', JSON.stringify(data,null,2));
  document.querySelector('#target').onchange = e => { target=e.target.value; stage=null; render(); };
  document.querySelectorAll('[data-tab]').forEach(b => b.onclick = () => {tab=b.dataset.tab; render();});
  document.querySelectorAll('[data-stage]').forEach(b => b.onclick = () => {stage=Number(b.dataset.stage); if(stage===1||stage===9)tab='evidence'; render();});
  document.querySelector('#close-stage')?.addEventListener('click',()=>{stage=null;render();});
  document.querySelector('#return-epsps')?.addEventListener('click',()=>{target='AT2G45300';render();});
}
function structureView() {
  const contexts = [...new Set(Object.values(structures).filter(s=>s.subject===subject).map(s=>s.context))];
  if(!contexts.includes(context))context=contexts[0];
  const model = Object.values(structures).find(s=>s.subject===subject&&s.context===context&&s.seed===seed);
  const replica = Object.values(structures).find(s=>s.subject===subject&&s.context===context&&s.seed!==seed);
  const labels = {herbicide:'Glyphosate + S3P',native:'PEP + S3P',native_s3p:'PEP + S3P · S3P score',s3p_only:'S3P only'};
  document.querySelector('#view').innerHTML = `<section class="structure-section"><div class="section-heading"><div><h2>Molecular structure</h2><span>${esc(subject==='arabidopsis_wt'?'Arabidopsis wild type':subject==='ecoli_wt'?'E. coli wild-type control':'E. coli G96A control')}</span></div><div class="viewer-selects"><label>Protein<select id="subject"><option value="arabidopsis_wt">Arabidopsis WT</option><option value="ecoli_wt">E. coli WT</option><option value="ecoli_G96A">E. coli G96A</option></select></label><label>Complex<select id="context">${contexts.map(c=>`<option value="${c}">${labels[c]}</option>`).join('')}</select></label><label>Seed<select id="seed">${[101,103,107].map(s=>`<option>${s}</option>`).join('')}</select></label></div></div>
      <div class="molecule-stage"><div id="molecule"></div><div class="scene-top"><span class="scene-label">BOLTZ-2 <span> / </span> ${model.protein.length} residues</span><div class="scene-actions">${iconButton('reset','scan','Fit full structure')}${iconButton('pocket','focus','Focus ligand pocket')}${iconButton('spin','rotate-3d','Rotate structure',spinning)}${iconButton('overlay','layers-2','Overlay another aligned seed',overlay)}</div></div><div class="scene-bottom"><div class="swatches"><span><i style="background:#218b77"></i>Protein trace</span><span><i style="background:#21718c"></i>Contact residues</span><span><i style="background:#d55249"></i>${context==='herbicide'?'Glyphosate':context==='s3p_only'?'S3P':'PEP'}</span>${context==='s3p_only'?'':'<span><i style="background:#d19a20"></i>S3P</span>'}<span id="overlay-key" ${overlay?'':'hidden'}><i style="background:#b95a78"></i>Seed ${replica.seed}</span></div><span id="residue-info">${Object.values(model.contacts).flat().filter((v,i,a)=>a.indexOf(v)===i).length} contact residues</span></div></div>
      ${legend('Figure 1 · Structure and ligand contacts','Actual predicted coordinates. The green tube follows protein C-alpha atoms; ligand spheres are heavy-atom positions, colored by ligand rather than element. Blue markers identify residues within 5 angstroms of any ligand heavy atom. Overlay is another seed of the same sequence and complex, globally aligned by C-alpha least squares, not an experimental reference. Arabidopsis numbering refers to the full sequence; E. coli uses observed-residue numbering. Shape agreement is not proof of retained activity.')}</section>
      <section class="results-section"><div class="section-heading"><div><h2>Current calibration panel</h2><span>Wild-type references and experimental homolog controls</span></div>${badge('Uncalibrated affinity predictions','warning')}</div><div class="table-scroll"><table><thead><tr><th>Protein / mutation</th><th>Glyphosate pIC50</th><th>PEP pIC50</th><th>S3P pIC50</th><th>Min TM-score</th><th>Max RMSD, Å</th></tr></thead><tbody>${data.summary.map(r=>`<tr><th>${esc(r.organism)} <b>${r.mutation}</b></th>${['glyphosate','pep','s3p'].map(l=>`<td>${r.ligands[l]?`${fmt(r.ligands[l].mean)}<small>${fmt(r.ligands[l].min)}–${fmt(r.ligands[l].max)}</small>`:'<span class="muted">Not scored</span>'}</td>`).join('')}<td>${fmt(r.tm_min)}</td><td>${fmt(r.rmsd_max)}</td></tr>`).join('')}</tbody></table></div>
      ${legend('Table 1 · Affinity predictions and structural checks','Values are mean predicted pIC50 and min–max across three seeds, not measured binding constants or confidence intervals. Higher pIC50 means stronger predicted binding. Minimum TM-score and maximum global C-alpha RMSD pool within-condition repeats and, for E. coli only, matched crystal comparisons; these are not mutant-to-WT structural similarity scores. E. coli G96A S3P is scored in S3P-only context, while Arabidopsis S3P is scored with PEP; they must not be compared directly. E. coli is a homolog control, not the Arabidopsis target.')}</section>`;
  for(const [id,value] of [['subject',subject],['context',context],['seed',seed]])document.querySelector('#'+id).value=value;
  document.querySelector('#subject').onchange=e=>{subject=e.target.value;render();};
  document.querySelector('#context').onchange=e=>{context=e.target.value;render();};
  document.querySelector('#seed').onchange=e=>{seed=Number(e.target.value);render();};
  try { viewer=createViewer(document.querySelector('#molecule'),model,replica,r=>{document.querySelector('#residue-info').textContent=`${r.name} ${r.position} · ligand contact`;}); }
  catch {document.querySelector('#molecule').innerHTML='<p class="viewer-error">3D rendering is unavailable in this browser. The verified structure results remain in the table below.</p>';return;}
  viewer.spin(spinning);viewer.overlay(overlay);
  document.querySelector('#reset').onclick=()=>viewer.reset();
  document.querySelector('#pocket').onclick=()=>viewer.pocket();
  document.querySelector('#spin').onclick=e=>{spinning=!spinning;viewer.spin(spinning);e.currentTarget.setAttribute('aria-pressed',spinning);};
  document.querySelector('#overlay').onclick=e=>{overlay=!overlay;viewer.overlay(overlay);e.currentTarget.setAttribute('aria-pressed',overlay);document.querySelector('#overlay-key').hidden=!overlay;};
}
function calibrationView() {
  const binding=data.binding;
  document.querySelector('#view').innerHTML=`<section class="calibration-summary"><div><span class="eyebrow">BINDING CALIBRATION</span><h2>${binding.status==='VALIDATED_WITHIN_SCOPE'?'Validated within assay scope':'Structure is repeatable. Binding is not yet calibrated.'}</h2><p>No experimental resistance or retained enzyme activity has been established.</p></div>${badge(pretty(binding.status),'warning')}</section>
    <section class="section"><div class="section-heading"><div><h2>Known-control check</h2><span>E. coli G96A compared with its wild-type reference</span></div></div><div class="control-grid">${binding.controls.map(c=>`<article class="control-item"><div class="control-title"><h3>${c.ligand==='pep'?'Native substrate · PEP':'Herbicide · Glyphosate'}</h3>${badge(c.status==='DIRECTION_CONSISTENT'?'Direction agrees':'Direction warning',c.status==='DIRECTION_CONSISTENT'?'success':'warning')}</div><div class="delta">${c.paired_delta_mean>0?'+':''}${fmt(c.paired_delta_mean)}<span>mutant − WT predicted pIC50</span></div><div class="direction-axis"><span>Weaker predicted</span><span>Stronger predicted</span><div class="axis-line"><i class="axis-zero"></i><i class="delta-marker" style="left:${Math.max(2,Math.min(98,50+c.paired_delta_mean*150))}%"></i></div></div><dl><div><dt>Wild-type mean</dt><dd>${fmt(c.wild_type_mean)}</dd></div><div><dt>G96A mean</dt><dd>${fmt(c.mutant_mean)}</dd></div><div><dt>Paired-seed difference range</dt><dd>${c.paired_delta_range.map(v=>fmt(v)).join(' to ')}</dd></div></dl><p class="control-note">${c.ligand==='pep'?'Published experiment: reduced PEP affinity. Model: stronger predicted binding.':'Published experiment: glyphosate insensitivity. Model: a small shift toward weaker predicted binding.'}</p><a href="${esc(c.source)}" target="_blank" rel="noreferrer">Experimental source ${icon('arrow-up-right')}</a></article>`).join('')}</div>${legend('Figure 2 · Directional binding-control diagnostic','Mean paired-seed differences use the same three seeds and matched homolog MSAs. Negative means weaker predicted binding. All three glyphosate differences are negative and all three PEP differences positive. Axis spans approximately -0.33 to +0.33 pIC50. Published G96A observations use kinetic endpoints, not model pIC50: these comparisons can flag direction disagreement but cannot calibrate a numerical scale. These E. coli controls are not independent Arabidopsis functional validation.',true)}</section>
    <section class="section"><div class="section-heading"><div><h2>Quantitative calibration readiness</h2><span>Endpoint-matched experiments are required before fitting</span></div><a href="/data/binding_calibration.json" download class="text-button">${icon('download')} Report</a></div><div class="calibration-steps">${[['01','Measured labels','IC50 only; no Kd, Ki or Km substitution'],['02','Independent splits','8 fit / 9 uncertainty / 10 test groups minimum'],['03','Held-out checks','Error, interval coverage and useful precision'],['04','Scope-limited output','No extrapolation or automatic function approval']].map(([n,title,body])=>`<div><span>${n}</span><h3>${title}</h3><p>${body}</p></div>`).join('')}</div><ul class="reasons">${binding.reasons.map(r=>`<li>${esc(r)}</li>`).join('')}</ul>${legend('Method legend · Calibration policy',binding.legend+' Minimum counts and acceptance limits are predeclared engineering safeguards, not a statistical power calculation or proof of biological validity.')}</section>
    <section class="section"><div class="section-heading"><div><h2>Structural calibration</h2><span>All 30 comparisons pass the original structural limits</span></div>${badge('30 / 30 checks','success')}</div><img class="scientific-figure" src="/data/structure-checks.png" alt="Measured structural deviations for predicted WT repeats and E. coli crystal controls, all below the 1 angstrom RMSD limit">${legend('Figure 3 · Structural matching',data.legends.structure)}</section>`;
}
function archiveView(){
  document.querySelector('#view').innerHTML=`<section class="section"><div class="section-heading"><div><h2>Earlier candidate screen</h2><span>Six point mutations + native WT reference · 19 September</span></div>${badge('Legacy protocol','warning')}</div><div class="notice">${icon('history')}<p>Query-only MSA, two seeds. All six candidates failed the computational screen, but WT variability was also high. These results do not isolate a mutation effect and are not directly comparable with the current matched-MSA run.</p></div><div class="table-scroll"><table><thead><tr><th>Mutation</th><th>Glyphosate pIC50</th><th>PEP pIC50</th><th>S3P pIC50</th><th>TM-score</th><th>Global RMSD, Å</th><th>Screen result</th></tr></thead><tbody>${data.archive.rows.map(r=>`<tr><th>${esc(r.mutation)}${r.mutation==='WT'?'<small>Native reference</small>':''}</th>${['glyphosate','pep','s3p'].map(l=>`<td>${fmt(r[l+'_pic50_mean'])}<small>${fmt(r[l+'_pic50_min'])}–${fmt(r[l+'_pic50_max'])}</small></td>`).join('')}<td>${fmt(r.query_tm_score)}</td><td>${fmt(r.global_ca_rmsd_angstrom)}</td><td>${badge(r.mutation==='WT'?'Reference':'Did not pass',r.mutation==='WT'?'neutral':'warning')}</td></tr>`).join('')}</tbody></table></div>${legend('Table 2 · Archived mutation screen',data.legends.archive,true)}</section>`;
}
function evidenceView(){
  const models=data.models;
  const judged=data.workflow_review?.judge?.status==='COMPLETED';
  document.querySelector('#view').innerHTML=`<section class="section"><div class="section-heading"><div><h2>Evidence and model access</h2><span>${esc(data.run.workflow_id || data.run.id)}</span></div>${badge(`${data.integrity.verified_artifacts} artifacts verified`,'success')}</div><div class="evidence-grid"><div><h3>Evidence model</h3><strong>${esc(models?.evidence?.model || 'Not configured')}</strong><p>${models?.evidence?.metadata_accessible || models?.evidence?.listed?'Model metadata accessible; synthesis status is reported separately.':'Account access unavailable. No substitute model used.'}</p>${badge('Synthesis: '+pretty(status('evidence')))}<p>Literature retrieval: ${esc(pretty(data.literature?.status || 'not_run'))}</p></div><div><h3>LLM judge</h3><strong>${esc(models?.judge?.model || 'Not configured')}</strong><p>${judged?'A real API request audited the workflow diagnostics and reviewer response. Candidate approval remains blocked.':models?.judge?.request_verified?'A standalone calibration review was verified; candidate judging remains separate.':'No live workflow judge response recorded.'}</p>${badge(judged?'Workflow judge completed':'Not verified',judged?'success':'warning')}</div></div>
    ${data.workflow_review?['review','judge'].map(role=>{const item=data.workflow_review[role];return `<div class="review"><div class="section-heading"><h3>Workflow ${role}</h3>${badge(pretty(item.status),item.status==='COMPLETED'?'success':'warning')}</div>${item.result?`<p>${esc(item.result.summary)}</p><h4>Next actions</h4><ul>${item.result.next_actions.map(x=>`<li>${esc(x)}</li>`).join('')}</ul>`:''}</div>`;}).join('')+legend('Diagnostic review legend',data.workflow_review.legend):''}
    ${models?.calibration_review?`<div class="review"><div class="section-heading"><h3>Standalone calibration review</h3>${badge(pretty(models.calibration_review.recommendation),'warning')}</div><h4>Uncertainty identified by the judge</h4><ul>${models.calibration_review.uncertainty.map(x=>`<li>${esc(x)}</li>`).join('')}</ul><p class="muted">This is an LLM assessment, not experimental evidence. It cannot override deterministic gates.</p></div>`:''}</section>
    <section class="section"><div class="section-heading"><h2>Source records</h2></div>${data.sources.map(s=>`<div class="source-row"><span class="source-id">${esc(s.id)}</span><div><strong>${esc(s.title || s.scope)}</strong><p>${s.evidence_level?'Retrieved abstract; applicability unreviewed':'Experimental crystal reference'} · source bytes hash-verified</p><code>${esc(s.sha256)}</code></div><a href="${esc(s.url)}" target="_blank" rel="noreferrer" aria-label="Open ${esc(s.id)} source">${icon('arrow-up-right')}</a></div>`).join('')}${legend('Evidence legend','1G6S is E. coli WT EPSPS with glyphosate and S3P. 1MI4 is the G96A homolog with S3P. These are homolog controls, not Arabidopsis experimental structures. '+(data.literature?.legend || 'No literature retrieval recorded.'))}</section>
    <section class="section"><div class="section-heading"><h2>Open scientific requirements</h2></div><ul class="requirements">${data.blocking_reasons.map(r=>`<li>${icon('circle-alert')}<span>${esc(r)}</span></li>`).join('')}</ul></section>`;
}
function render(){viewer?.dispose();viewer=null;shell();if(target==='AT2G45300')({structure:structureView,calibration:calibrationView,archive:archiveView,evidence:evidenceView}[tab])();createIcons({icons});}
async function start(){
  try{
    const responses=await Promise.all(['/data/results.json','/data/structures.json'].map(url=>fetch(url)));
    if(responses.some(r=>!r.ok))throw new Error('Missing exported data');
    [data,structures]=await Promise.all(responses.map(r=>r.json()));
    render();
  }catch(error){root.innerHTML=`<main class="empty"><h1>Scientific results unavailable</h1><p>The recorded-run data could not be loaded. No substitute results have been generated.</p><button id="retry" class="primary">Retry</button></main>`;document.querySelector('#retry').onclick=start;}
}
start();
