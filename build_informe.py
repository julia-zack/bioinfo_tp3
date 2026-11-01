"""informe.md -> informe.html -> informe.pdf, con el mismo formato que los
informes de los TP anteriores: una sola pagina HTML con las figuras adentro
(en base64) y el pdf impreso desde Chrome sin interfaz.

    python3 build_informe.py

Necesita `pip install markdown` y Google Chrome instalado.
"""

import base64
import os
import re
import subprocess

import markdown

HERE = os.path.dirname(os.path.abspath(__file__))
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

CSS = """
:root {
  --ink: #1a1a1a;
  --ink-soft: #4a4a4a;
  --rule: #e0ddd7;
  --paper: #fdfdfc;
  --code-bg: #f4f2ee;
  --accent: #2e5c8a;
}

html { -webkit-text-size-adjust: 100%; }

body {
  margin: 0 auto;
  padding: 4rem 1.5rem 6rem;
  max-width: 46rem;
  background: var(--paper);
  color: var(--ink);
  font: 16px/1.7 Charter, "Iowan Old Style", Georgia, "Times New Roman", serif;
}

h1 {
  font-size: 2rem;
  line-height: 1.25;
  margin: 0 0 2.5rem;
  padding-bottom: 1rem;
  border-bottom: 2px solid var(--ink);
  font-weight: 600;
}

h2 {
  font-size: 1.5rem;
  margin: 3.5rem 0 1.25rem;
  padding-bottom: 0.4rem;
  border-bottom: 1px solid var(--rule);
  font-weight: 600;
}

h3 {
  font-size: 1.15rem;
  margin: 2.5rem 0 0.9rem;
  font-weight: 600;
  color: var(--ink-soft);
}

h4 {
  font-size: 1rem;
  margin: 2rem 0 0.8rem;
  font-weight: 600;
  color: var(--ink-soft);
}

p { margin: 0 0 1.15rem; }

ul, ol { margin: 0 0 1.15rem; padding-left: 1.4rem; }
li { margin-bottom: 0.4rem; }

blockquote {
  margin: 0 0 1.75rem;
  padding: 1rem 1.25rem;
  border-left: 3px solid var(--accent);
  background: #f7f9fb;
  color: var(--ink-soft);
  font-size: 0.95rem;
}
blockquote p:last-child { margin-bottom: 0; }

pre {
  margin: 0 0 1.75rem;
  padding: 1rem 1.15rem;
  white-space: pre-wrap;
  overflow-wrap: anywhere;
  background: var(--code-bg);
  border: 1px solid var(--rule);
  border-radius: 4px;
  font: 12px/1.5 "SF Mono", Menlo, Consolas, monospace;
}
pre code { background: none; padding: 0; font-size: inherit; }

code {
  padding: 0.12em 0.35em;
  background: var(--code-bg);
  border-radius: 3px;
  font: 0.88em/1 "SF Mono", Menlo, Consolas, monospace;
}

table {
  width: 100%;
  margin: 0 0 1.75rem;
  border-collapse: collapse;
  font-size: 0.92rem;
}
th, td {
  padding: 0.5rem 0.7rem;
  border-bottom: 1px solid var(--rule);
  text-align: left;
}
th {
  border-bottom: 2px solid var(--ink-soft);
  font-weight: 600;
}
td:not(:first-child), th:not(:first-child) { text-align: right; }

figure {
  margin: 2rem 0;
  text-align: center;
}
figure img {
  max-width: 100%;
  height: auto;
  border: 1px solid var(--rule);
  border-radius: 4px;
}
figcaption {
  margin-top: 0.6rem;
  color: var(--ink-soft);
  font-size: 0.85rem;
  font-style: italic;
}

hr {
  margin: 3rem 0;
  border: 0;
  border-top: 1px solid var(--rule);
}

a { color: var(--accent); }
strong { font-weight: 600; }
em { color: var(--ink-soft); }

@media print {
  body { padding: 0; max-width: none; font-size: 10.5pt; }
  h1 { margin-bottom: 1.25rem; }
  h2 { margin-top: 2rem; }
  h2, h3, h4 { break-after: avoid; page-break-after: avoid; }
  blockquote, figure, table {
    break-inside: avoid; page-break-inside: avoid;
  }
  th, td { padding: 0.2rem 0.6rem; }
  table { margin-bottom: 1.2rem; font-size: 0.9em; }
  hr { display: none; }
  .page-break { break-before: page; page-break-before: always; }
}
"""


def figure(match):
    """![pie](ruta) -> <figure> con la imagen en base64 y el pie como caption."""
    caption, path = match.group(1), match.group(2)
    with open(os.path.join(HERE, path), "rb") as f:
        data = base64.b64encode(f.read()).decode()
    return (f'<figure><img alt="{caption}" src="data:image/png;base64,{data}">'
            f'<figcaption>{caption}</figcaption></figure>')


def main():
    with open(os.path.join(HERE, "informe.md")) as f:
        text = f.read()
    title = re.match(r"# (.*)", text).group(1)

    body = markdown.markdown(text, extensions=["tables", "fenced_code"])
    # python-markdown envuelve cada imagen en un <p>; <figure> no va dentro de un <p>
    body = re.sub(r'<p><img alt="([^"]*)" src="([^"]*)" ?/?></p>', figure, body)

    html = (f'<!doctype html>\n<html lang="es">\n<head>\n<meta charset="utf-8">\n'
            f'<meta name="viewport" content="width=device-width, initial-scale=1">\n'
            f'<title>{title}</title>\n<style>{CSS}</style>\n</head>\n<body>\n{body}\n</body>\n</html>\n')
    out_html = os.path.join(HERE, "informe.html")
    with open(out_html, "w") as f:
        f.write(html)

    out_pdf = os.path.join(HERE, "informe.pdf")
    subprocess.run([CHROME, "--headless", "--disable-gpu", "--no-pdf-header-footer",
                    f"--print-to-pdf={out_pdf}", f"file://{out_html}"],
                   check=True, capture_output=True)
    print(f"{out_html}\n{out_pdf}")


if __name__ == "__main__":
    main()
