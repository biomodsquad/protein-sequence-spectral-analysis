"""Check built HTML's local links, asset paths, and fragment targets without network."""

import argparse
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


class Page(HTMLParser):
    """Collect render-relevant links and anchors from one generated page."""

    def __init__(self, text):
        super().__init__()
        self.links, self.anchors = [], set()
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        """Record assets/navigation and named anchors; external URLs are filtered later."""
        attrs = dict(attrs)
        if 'id' in attrs:
            self.anchors.add(attrs['id'])
        if tag == 'a' and 'name' in attrs:
            self.anchors.add(attrs['name'])
        attribute = {'a': 'href', 'link': 'href', 'img': 'src', 'script': 'src',
                     'form': 'action'}.get(tag)
        if attribute and attrs.get(attribute):
            self.links.append(attrs[attribute])


def main():
    """Fail on broken local links; external network links are deliberately not checked."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    root = parser.parse_args().directory.resolve()
    pages = {p.resolve(): Page(p.read_text()) for p in root.rglob('*.html')}
    if not pages:
        raise ValueError('No built HTML pages found')
    broken, checked = [], 0
    for source, page in pages.items():
        for target in page.links:
            url = urlsplit(target)
            if url.scheme or url.netloc:
                continue
            path = unquote(url.path)
            resolved = ((root / path.lstrip('/')) if path.startswith('/')
                        else source.parent / path).resolve() if path else source
            if resolved.is_dir():
                resolved = resolved / 'index.html'
            checked += 1
            if not resolved.is_file():
                broken.append(f'{source.relative_to(root)}: missing {target}')
            elif url.fragment and resolved in pages and unquote(url.fragment) not in pages[resolved].anchors:
                broken.append(f'{source.relative_to(root)}: missing anchor {target}')
    if broken:
        raise RuntimeError('\n'.join(sorted(set(broken))))
    print(f'Checked {checked} local links/assets across {len(pages)} HTML pages')


if __name__ == '__main__':
    main()
