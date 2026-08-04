// ---- Newick parser ----
function parseNewick(str) {
  str = str.trim().replace(/;\s*$/, '');
  let pos = 0;
 
  function parseNode() {
    const node = { children: [] };
    if (str[pos] === '(') {
      pos++; // skip (
      node.children.push(parseNode());
      while (str[pos] === ',') {
        pos++;
        node.children.push(parseNode());
      }
      pos++; // skip )
    }
    let name = '';
    while (pos < str.length && !',()'.includes(str[pos]) && str[pos] !== ':') {
      name += str[pos++];
    }
    node.name = name;
    if (str[pos] === ':') {
      pos++;
      let len = '';
      while (pos < str.length && !',()'.includes(str[pos])) {
        len += str[pos++];
      }
      node.length = parseFloat(len);
    }
    return node;
  }
 
  return parseNode();
}
 
// ---- Layout ----
function layoutTree(root) {
  let leafIndex = 0;
  function assignY(node) {
    if (!node.children.length) {
      node.y = leafIndex++;
    } else {
      node.children.forEach(assignY);
      const ys = node.children.map(c => c.y);
      node.y = (Math.min(...ys) + Math.max(...ys)) / 2;
    }
  }
  assignY(root);
 
  function assignX(node, parentX) {
    const len = (node.length !== undefined && !isNaN(node.length)) ? node.length : 1;
    node.x = parentX + len;
    node.children.forEach(c => assignX(c, node.x));
  }
  root.x = 0;
  root.children.forEach(c => assignX(c, 0));
 
  return leafIndex; // total leaf count
}
 
// ---- Render ----
function renderTree(newickStr, svgId) {
  const root = parseNewick(newickStr);
  const leafCount = layoutTree(root);
 
  const leafSpacing = 40;
  const marginLeft = 20;
  const marginRight = 160; // room for leaf labels
  const marginTop = 30;
  const marginBottom = 30;
 
  let maxX = 0;
  (function findMaxX(n) { maxX = Math.max(maxX, n.x); n.children.forEach(findMaxX); })(root);
 
  const xScale = maxX > 0 ? (600 / maxX) : 1;
  const width = marginLeft + 600 + marginRight;
  const height = marginTop + leafCount * leafSpacing + marginBottom;
 
  const svg = document.getElementById(svgId);
  svg.setAttribute('width', width);
  svg.setAttribute('height', height);
  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  svg.innerHTML = '';
 
  function px(node) { return marginLeft + node.x * xScale; }
  function py(node) { return marginTop + node.y * leafSpacing; }
 
  function draw(node, parent) {
    if (parent) {
      // horizontal branch from parent.x to node.x, at node.y
      const line1 = document.createElementNS('http://www.w3.org/2000/svg', 'path');
      line1.setAttribute('class', 'branch');
      line1.setAttribute('d', `M ${px(parent)} ${py(node)} H ${px(node)}`);
      svg.appendChild(line1);
    }
 
    if (node.children.length) {
      // vertical connector spanning children's y at node.x
      const ys = node.children.map(py);
      const line2 = document.createElementNS('http://www.w3.org/2000/svg', 'path');
      line2.setAttribute('class', 'branch');
      line2.setAttribute('d', `M ${px(node)} ${Math.min(...ys)} V ${Math.max(...ys)}`);
      svg.appendChild(line2);
 
      // internal node label + dot
      if (node.name) {
        const dot = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
        dot.setAttribute('class', 'node-dot');
        dot.setAttribute('cx', px(node));
        dot.setAttribute('cy', py(node));
        dot.setAttribute('r', 3);
        svg.appendChild(dot);
 
        const lx = px(node) + 6;
        const ly = py(node) + 8;
        const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        label.setAttribute('class', 'node-label');
        label.setAttribute('x', lx);
        label.setAttribute('y', ly);
        label.setAttribute('transform', `rotate(+45 ${lx} ${ly})`);
        label.textContent = node.name;
        svg.appendChild(label);
      }
 
      node.children.forEach(c => draw(c, node));
    } else {
      // leaf label
      const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
      label.setAttribute('class', 'leaf-label');
      label.setAttribute('x', px(node) + 8);
      label.setAttribute('y', py(node));
      label.textContent = node.name;
      svg.appendChild(label);
    }
  }
 
  draw(root, null);
}
 
// ---- Wire up one tree instance (textarea + svg pair) ----
function initTree(textareaId, svgId) {
  const textarea = document.getElementById(textareaId);
  renderTree(textarea.value, svgId);
  textarea.addEventListener('input', () => {
    try {
      renderTree(textarea.value, svgId);
    } catch (e) {
      // ignore parse errors while user is mid-edit
    }
  });
}
