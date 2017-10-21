import urllib.request


def connect_url(url, decode):
    """
    :method: a method to perfect urllib.request.connect()
    :param url: string
    :param decode: string like "utf-8"
    :return: string(the url content)
    """
    req = urllib.request.Request(url)
    response = urllib.request.urlopen(req)
    html = response.read()
    html = html.decode(decode)
    return html


def callback(block, block_size, total_size):
    """
    :method: a method to print download schedule
    :param block:
    :param block_size:
    :param total_size:
    :return:
    """
    percent = 100.0 * block * block_size / total_size
    if percent > 100:
        percent = 100
    print("%.2f%%" % percent)
