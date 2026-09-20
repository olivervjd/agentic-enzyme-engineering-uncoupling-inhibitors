import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

export function createViewer(container, model, replica, onResidue) {
  const scene = new THREE.Scene();
  scene.background = new THREE.Color('#f1f6f4');
  const renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  container.append(renderer.domElement);
  renderer.domElement.setAttribute('aria-label', `${model.subject}, ${model.context}, seed ${model.seed}: interactive molecular structure`);
  renderer.domElement.setAttribute('role', 'img');
  renderer.domElement.dataset.residues = model.protein.length;
  const camera = new THREE.PerspectiveCamera(38, 1, .1, 2000);
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.autoRotateSpeed = .65;
  controls.minDistance = 4;
  controls.maxDistance = 220;
  scene.add(new THREE.HemisphereLight(0xffffff, 0x7c9b8e, 2.5));
  const light = new THREE.DirectionalLight(0xffffff, 3);
  light.position.set(30, 60, 80);
  scene.add(light);
  const group = new THREE.Group();
  scene.add(group);
  const overlay = new THREE.Group();
  group.add(overlay);
  overlay.visible = false;
  const points = model.protein.map(r => new THREE.Vector3(...r.xyz));
  const bounds = new THREE.Box3().setFromPoints(points);
  const center = bounds.getCenter(new THREE.Vector3());
  const radius = bounds.getSize(new THREE.Vector3()).length() / 2;
  const clickable = [];
  const materials = [];
  const material = (color, opacity=1) => {
    const value = new THREE.MeshStandardMaterial({ color, roughness: .42, metalness: .08, transparent: opacity < 1, opacity });
    materials.push(value);
    return value;
  };
  function trace(rows, color, parent, thickness=.42, opacity=1) {
    const mat = material(color, opacity);
    let segment = [];
    const flush = () => {
      if (segment.length > 1) {
        const curve = new THREE.CatmullRomCurve3(segment);
        parent.add(new THREE.Mesh(new THREE.TubeGeometry(curve, segment.length * 5, thickness, 7, false), mat));
      }
      segment = [];
    };
    for (const row of rows) {
      const point = new THREE.Vector3(...row.xyz);
      if (segment.length && segment.at(-1).distanceTo(point) > 6) flush();
      segment.push(point);
    }
    flush();
  }
  trace(model.protein, '#218b77', group);
  if (replica) trace(replica.protein, '#b95a78', overlay, .28, .72);
  const contactPositions = new Set(Object.values(model.contacts).flat());
  const contactMaterial = material('#21718c');
  const sphere = new THREE.SphereGeometry(1, 12, 10);
  for (const row of model.protein) {
    if (!contactPositions.has(row.position)) continue;
    const mesh = new THREE.Mesh(sphere, contactMaterial);
    mesh.position.set(...row.xyz);
    mesh.scale.setScalar(.85);
    mesh.userData = row;
    group.add(mesh);
    clickable.push(mesh);
  }
  const ligandColors = { B: '#d55249', C: '#d19a20' };
  for (const ligand of model.ligands) {
    const mat = material(ligandColors[ligand.chain] || '#905eac');
    for (const atom of ligand.atoms) {
      const mesh = new THREE.Mesh(sphere, mat);
      mesh.position.set(...atom.xyz);
      mesh.scale.setScalar(atom.element === 'P' ? 1.15 : .8);
      group.add(mesh);
    }
  }
  function focus(target, size) {
    controls.target.copy(target);
    const fov = THREE.MathUtils.degToRad(camera.fov);
    const distance = size / Math.sin(Math.atan(Math.tan(fov / 2) * Math.min(camera.aspect, 1))) * 1.08;
    camera.position.copy(target).add(new THREE.Vector3(.6, .25, 1).normalize().multiplyScalar(distance));
    controls.update();
  }
  const resize = new ResizeObserver(() => {
    const {width, height} = container.getBoundingClientRect();
    if (!width || !height) return;
    renderer.setSize(width, height);
    camera.aspect = width / height;
    camera.updateProjectionMatrix();
    focus(center, radius);
  });
  resize.observe(container);
  const raycaster = new THREE.Raycaster();
  const pointer = new THREE.Vector2();
  const click = event => {
    const rect = renderer.domElement.getBoundingClientRect();
    pointer.set((event.clientX - rect.left) / rect.width * 2 - 1, -(event.clientY - rect.top) / rect.height * 2 + 1);
    raycaster.setFromCamera(pointer, camera);
    const hit = raycaster.intersectObjects(clickable)[0];
    if (hit) onResidue(hit.object.userData);
  };
  renderer.domElement.addEventListener('click', click);
  let frame;
  function animate() { frame = requestAnimationFrame(animate); controls.update(); renderer.render(scene, camera); }
  animate();
  return {
    reset: () => focus(center, radius),
    pocket: () => {
      const atoms = model.ligands.flatMap(l => l.atoms.map(a => new THREE.Vector3(...a.xyz)));
      if (!atoms.length) return;
      const box = new THREE.Box3().setFromPoints(atoms);
      focus(box.getCenter(new THREE.Vector3()), Math.max(10, box.getSize(new THREE.Vector3()).length() / 2 + 4));
    },
    spin: value => { controls.autoRotate = value; },
    overlay: value => { overlay.visible = value; },
    dispose: () => {
      cancelAnimationFrame(frame); resize.disconnect(); controls.dispose();
      renderer.domElement.removeEventListener('click', click);
      scene.traverse(item => { if (item.geometry) item.geometry.dispose(); });
      materials.forEach(m => m.dispose()); renderer.dispose(); renderer.domElement.remove();
    }
  };
}
