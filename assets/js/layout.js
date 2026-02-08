document.addEventListener("DOMContentLoaded", () => {
  const isInPages = window.location.pathname.includes("/pages/");
  const basePath = isInPages ? ".." : ".";
  const target = document.getElementById("site-header");

  if (!target) return;

  fetch(`${basePath}/assets/partials/header.html`)
    .then((res) => res.text())
    .then((html) => {
      const rendered = html.replace(/{{BASE}}/g, basePath);
      target.innerHTML = rendered;

      const path = window.location.pathname;
      const navMap = {
        "/pages/about.html": "about",
        "/pages/guide.html": "guide",
        "/pages/institute.html": "institute",
        "/pages/publications.html": "publications",
      };
      const activeKey = navMap[path];

      const links = target.querySelectorAll(".navlinks a");
      links.forEach((link) => {
        if (activeKey && link.dataset.nav === activeKey) {
          link.classList.add("active");
        } else {
          link.classList.remove("active");
        }
      });
    })
    .catch(() => {});
});
