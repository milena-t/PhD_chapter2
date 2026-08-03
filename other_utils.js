function openTab(evt, cityName) {
    // only change tabs on the same level to allow nested tabs etc.
    const tab = document.getElementById(cityName);
    const parenttab = tab.parentNode;
  
    // Hide sibling tab panels.
    for (const el of document.getElementsByClassName("tabcontent")) {
      if (el.parentNode === parenttab) el.style.display = "none";
    }
  
    // Clear "active" only on tablinks at the same level. A tablink's grandparent
    // (button -> .tab div -> parenttab) is the panel container.
    for (const btn of document.getElementsByClassName("tablinks")) {
      if (btn.parentNode.parentNode === parenttab) btn.classList.remove("active");
    }
  
    tab.style.display = "block";
    evt.currentTarget.classList.add("active");
}
  
document.addEventListener("DOMContentLoaded", function () {
    document.querySelectorAll(".tab").forEach(function (tabBar) {
      var firstBtn = tabBar.querySelector("button.tablinks");
      if (firstBtn) firstBtn.click();
    });
});

function loadCSVTable(csvPath, elementId) {
  Papa.parse(csvPath, {
    download: true,
    header: false,
    complete: function(results) {
      const rows = results.data;
      let html = "<table style='width:100%; border-collapse:collapse;'>";

      rows.forEach(row => {
        html += "<tr>";
        row.forEach(cell => {
          html += `<td style="border:1px solid #ccc; padding:4px;">${cell}</td>`;
        });
        html += "</tr>";
      });

      html += "</table>";
      document.getElementById(elementId).innerHTML = html;
    }
  });
}

