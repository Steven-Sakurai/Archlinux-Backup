import os
import requests
import re
from sys import argv
script, tag, begin, end = argv

import time
start_time = time.time()


# download single picture
# url: https://.../.../.../...jpg
def download(folder, url):
    if not os.path.exists(folder):
        os.makedirs(folder)
    req = requests.get(url)
    if req.status_code == requests.codes.ok:
        name = url.split('/')[-1]  
        f = open("./"+folder+'/'+name,'wb')
        f.write(req.content)
        f.close()
        return True
    else:
        return False

header = {'User-Agent':'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_12_5) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/60.0.3112.72 Safari/537.36'}
# url here: the contents page
# ex.: 'https://yande.re/post?tags=seifuku'
def fetch(url):
    r = requests.get(url,headers=header)
    text= r.text
    imgs=[]
    jpeg = re.compile(r'https://files.yande.re/image/[^\s]*?\.jpg')
    png = re.compile(r'https://files.yande.re/image/[^\s]*?\.png')

    imgs += jpeg.findall(text)
    imgs += png.findall(text)
    errors = []

    folder = url.split('=')[-1]
    folder = 'yandre-' + folder
    for img_url in imgs:
        if download(folder,img_url):
            print("download :"+img_url)
        else:
            errors.append(img_url)
    return errors



errs = []
if int(begin)==1:
    begin = int(begin)+1
    url1 = 'https://yande.re/post?tags=' + tag
    print(url1)
    errs += fetch(url1)
    print("------")
    print("Download images on %s complete~"%(url1))
    print("------")
    print("--- Total Time elapsed: %s seconds ---" % (time.time() - start_time))
    print("------")
if int(end)>1:
    for i in range(int(begin), int(end)+1):
        url = 'https://yande.re/post?page=' + str(i) + '&tags=' + tag
        print(url)
        errs += fetch(url)
        print("------")
        print("Download images on %s complete~"%(url))
        print("------")
        print("--- Total Time elapsed: %s seconds ---" % (time.time() - start_time))
        print("------")

print("All Download Tasks Complete~~~")
print("ERROR URLS:")
print(errs)







