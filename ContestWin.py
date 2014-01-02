# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from selenium import webdriver
from selenium.common.exceptions import TimeoutException, StaleElementReferenceException
from selenium.webdriver.support.ui import WebDriverWait 
import bs4
import urllib2
import time

# <codecell>

def get_current_artist():
    station_url = 'http://p2.wrff-fm.ccomrcdn.com/player/player_dispatcher.html?section=radio&action=listen_live'
    data = urllib2.urlopen(station_url).read()
    soup = bs4.BeautifulSoup(data)
    return soup.playercontent.justplayed.song.artist.attrs['name']
print get_current_artist()

# <codecell>

def send_txt_message(chrome_driver, txt_number, message):
    new_message_button = chrome_driver.find_element_by_xpath('//*[@id="newSms"]/div[2]')
    new_message_button.click()
    
    to_field = driver.find_element_by_xpath('//*[@id="selectContactForSingleCompose"]')
    to_field.send_keys(str(txt_number))
    
    mes_field = driver.find_element_by_xpath('//*[@id="send-one-text"]')
    mes_field.click()
    mes_field.send_keys(message)
    
    send_button = driver.find_element_by_xpath('//*[@id="send-button-single-text"]')
    send_button.click()

# <codecell>

#profile = webdriver.firefox.firefox_profile.FirefoxProfile('/home/will/.mozilla/firefox/fsg0yfdg.default/')
driver = webdriver.Chrome()

# <codecell>

driver.get('http://mightytext.net')

# <codecell>

tmp = driver.find_element_by_link_text('Login')
tmp.click()
time.sleep(10)

# <codecell>

ntmp = driver.find_element_by_xpath('//*[@id="Email"]')
ntmp.send_keys('judowill')
time.sleep(10)
otmp = driver.find_element_by_xpath('//*[@id="Passwd"]')
otmp.send_keys('judo8675309')
time.sleep(10)
rtmp = driver.find_element_by_xpath('//*[@id="signIn"]')
rtmp.click()
time.sleep(10)

# <codecell>


last_artist = ''
warned = False
while True:
    time.sleep(1)
    artist = get_current_artist()
    if artist != last_artist:
        print artist
        last_artist = artist
    if 'muse' in artist.lower():
        print 'PLAYING SONG!'
        send_txt_message(driver, 91045, 'Muse')
        time.sleep(5)
        if not warned:
            print 'telling cat!'
            try:
                send_txt_message(driver, 2157405170, 'playing the Muse, send messages now!!!')
            except:
                pass
            warned = True
            time.sleep(5)
        
        time.sleep(60*10)
        warned = False

# <codecell>


# <codecell>


