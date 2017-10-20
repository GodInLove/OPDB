import urllib.request


def connect_url(url, decode):
    req = urllib.request.Request(url)
    response = urllib.request.urlopen(req)
    html = response.read()
    html = html.decode(decode)
    return html


def callback(block, block_size, total_size):
    percent = 100.0 * block * block_size / total_size
    if percent > 100:
        percent = 100
    print("%.2f%%" % percent)
