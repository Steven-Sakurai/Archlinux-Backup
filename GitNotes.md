### Speed up git in China with SSR

```bash
git config --global http.postBuffer 2000000000
git config --global core.compression -1
git config --global http.proxy socks5://127.0.0.1:1080 
git config --global https.proxy socks5://127.0.0.1:1080 
```
