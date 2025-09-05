/*!
This is a patched 'Cookie helper functions' from Doxygen that overrides all
cookie setting and places everything in sessionStorage. It might break some browsers but works in Firefox and Chrome.

Original header:
 Copyright (c) 2023 Dimitri van Heesch
 Released under MIT license.
*/
let Cookie = {
  cookie_namespace: 'doxygen_',

  readSetting(cookie,defVal) {
    const val = sessionStorage.getItem(this.cookie_namespace+cookie);
    if (val) return val;
    return defVal;
  },

  writeSetting(cookie,val,days=0) { // default days='session', 0=session cookie, -1=delete
    sessionStorage.setItem(this.cookie_namespace+cookie,val);
  },

  eraseSetting(cookie) {
    sessionStorage.removeItem(this.cookie_namespace+cookie);
  },
}
