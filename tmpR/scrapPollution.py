
# Packages for scrapping
from selenium import webdriver
from selenium.webdriver.common.by import By
import pandas as pd
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import urllib.request

import time
import random


def main():
    # Selenium Options
    options = Options()
    #   enable headless mode
    # options.add_argument('--headless=new')
    # options.add_argument('window-size=1920x1080')
    #options.add_argument("--start-maximized")
    driver = webdriver.Edge()#executable_path=r'C:\Users\j53vande\Desktop\Webscrapping\msedgedriver.exe')
    driver.maximize_window()
    
    # General Vars
    EXPLICIT_WAIT_TIME = 10
    data = pd.DataFrame(
        columns=['date','12a', '1a', '2a','3a', '4a', '5a', '6a', 
                 '7a', '8a','9a', '10a', '11a', '12p', '1p', 
                 '2p', '3p', '4p', '5p', '6p', 
                 '7p', '8p', '9p', '10p', '11p'])
    
    
    driver.get('https://www.airqualityontario.com/history/pollutant.php?stationid=31129&pol_code=36&start_day=1&start_month=1&start_year=2007&showType=table&station_id=31129')
    
    while True:
        time.sleep(2*random.random())
        WebDriverWait(driver, EXPLICIT_WAIT_TIME).until(
            EC.frame_to_be_available_and_switch_to_it((By.XPATH,"//iframe[@id='chart-frame']"))
        )

        table = driver.find_element(By.CLASS_NAME,'resourceTable')
        rows = table.find_elements(By.XPATH,'./tbody/tr')
        
        for r in rows:
            data.loc[len(data)] = pd.Series(['']*25)
            
            elements = r.find_elements(By.XPATH,'./td')
            for e_idx in range(len(elements)):
                data.iloc[len(data)-1,e_idx] = elements[e_idx].text
                
        try:
            url = driver.find_element(By.LINK_TEXT,'Next').get_attribute('href')
            driver.switch_to.default_content()
            driver.get(url)
        except:
            break
    
    data.to_csv('C:/Users/jerem/Downloads/nox.csv')
        
        
    
if __name__ == '__main__':
    res = main()
    