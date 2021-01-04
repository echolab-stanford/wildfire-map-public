source("work/00_functions.R")


########################################################################################
# Written by: Anne Driscoll
# Downloads airport data
########################################################################################

months = c("January", "February", "March", "April", "May", "June", "July", "August", 
           "September", "October", "November", "December")
url = "https://www.transtats.bts.gov/Tables.asp?DB_ID=120&DB_Name=Airline%20On-Time%20Performance%20Data&DB_Short_Name=On-Time"

rD = rsDriver(geckover="0.23.0")
remDr = rD[["client"]]
remDr$navigate(url)
t = max(abs(rnorm(1, mean=3.5, sd=1.5)), 2)
Sys.sleep(t)

option = remDr$findElement(using = 'css selector', value = "tr:nth-child(7) a~ a+ a")
option$clickElement()

boxes = c("Year", "Month", "OriginAirportSeqID", "DestAirportSeqID", 
          "OriginAirportID", "DestAirportID", "Origin", "Dest", "TaxiOut", "TaxiIn")
for (box in boxes) {
    option = remDr$findElement(using='xpath', value=paste0('//input[@title="', box, '"]'))
    option$clickElement()
    Sys.sleep(abs(rnorm(1, mean=.6, sd=0.2)))
}


for (year in years) {
    for (month in months) {
        cat("Working on year:", year, ",", month, "\n")
        
        option = remDr$findElement(using = 'xpath', value = paste0("//*/option[text()=", year, "]"))
        option$clickElement()
        
        option = remDr$findElement(using = 'xpath', value = paste0('//*/option[text()="', month, '"]'))
        option$clickElement()
        
        Sys.sleep(abs(rnorm(1, mean=.2, sd=0.1)))
        
        download = remDr$findElement(using='xpath', value='//button[@name="Download"]')
        download$clickElement()
        
        while (!any(grepl("_T_.+zip$", list.files("~/Downloads")))) { #while downloading
            Sys.sleep(1)
        }
        
        unzip_rename(paste0(month, "_", year, "_ontime"))
        
        Sys.sleep(abs(rnorm(1, mean=.5, sd=0.1)))
    }
}


url = "https://www.transtats.bts.gov/Tables.asp?DB_ID=125&DB_Name=Airline%20Origin%20and%20Destination%20Survey%20%28DB1B%29&DB_Short_Name=Origin%20and%20Destination%20Survey"
remDr$navigate(url)
Sys.sleep(t)

option = remDr$findElement(using = 'css selector', value = "tr:nth-child(10) a~ a+ a")
option$clickElement()

boxes = c("Year", "Quarter", "OriginAirportSeqID", "OriginCityMarketID", "Origin", "Passengers")
for (box in boxes) {
    option = remDr$findElement(using='xpath', value=paste0('//input[@title="', box, '"]'))
    option$clickElement()
    Sys.sleep(abs(rnorm(1, mean=.6, sd=0.2)))
}

for (year in years) {
    for (q in c("Q 1", "Q 2", "Q 3", "Q 4")) {
        cat("Working on year:", year, ",", q, "\n")
        
        option = remDr$findElement(using = 'xpath', value = paste0("//*/option[text()=", year, "]"))
        option$clickElement()
        
        option = remDr$findElement(using = 'xpath', value = paste0('//*/option[text()="', q, '"]'))
        option$clickElement()
        
        Sys.sleep(abs(rnorm(1, mean=.2, sd=0.1)))
        
        download = remDr$findElement(using='xpath', value='//button[@name="Download"]')
        download$clickElement()
        
        while (!any(grepl("_T_.+zip$", list.files("~/Downloads")))) { #while downloading
            Sys.sleep(1)
        }
        
        x = unzip_rename(paste0(q, "_", year, "_tickets"))
        print(x)
        
        Sys.sleep(abs(rnorm(1, mean=1, sd=0.5)))
    }
}
