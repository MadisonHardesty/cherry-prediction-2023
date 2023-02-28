library("tidyverse")
library("rnoaa")
library("lubridate")
library("TDPanalysis")
library("caret")
library("glmnet")
library("forecast")

# import data
cherry <- read.csv("data/washingtondc.csv") %>% 
  bind_rows(read.csv("data/liestal.csv")) %>% 
  bind_rows(read.csv("data/kyoto.csv"))
cherry <- cherry[-216,] # remove kyoto row for 2023
# View(cherry)

# according to recommendations by Anderson et al. (1986)
# I kept finding base temperatures between 3 and 8 so I used the avergage for this model
base_temp <- 5.5 # under this temp, development stops. 

# get temperature data & calculate GDD based on the average method and Baskerville-Ermin method
stations <- ghcnd_stations()

get_temperature <- function (stationid, lat) {
  ghcnd_search(stationid = stationid, var = c("tmax", "tmin"), 
               date_min = "1900-01-01", date_max = "2022-12-31") %>%
    reduce(left_join) %>%
    transmute(year = parse_number(format(date, "%Y")),
              date, doy = date.to.DOY(format(date, "%Y/%m/%d"), format = 'yyyy/mm/dd'),
              tmax = tmax / 10, tmin = tmin / 10,
              avg_temp = (tmax + tmin) / 2,
              W = (tmax - tmin)/2,
              A = asin((base_temp - avg_temp) / W),
              gdd = case_when(
                ((tmax < base_temp) | (is.na(tmin) == TRUE) | (is.na(tmax) == TRUE)) ~ 0,
                (tmin >= base_temp) ~ avg_temp - base_temp,
                (tmin < base_temp) ~ ((W * cos(A)) - ((base_temp - avg_temp) * ((3.14/2)-A)))/3.14),
              cd = case_when(
                (tmax < 7.2) ~ 1,
                .default = 0)
    )
}

historic_temperatures <-
  tibble(location = "washingtondc", get_temperature("USC00186350", 38.88535)) %>%
  bind_rows(tibble(location = "liestal", get_temperature("GME00127786", 47.48140))) %>%
  bind_rows(tibble(location = "kyoto", get_temperature("JA000047759", 35.01198))) 
# View(historic_temperatures)

# calculate cumulative sum of GDD (growing degree days) and GDD^2 and GDD^3 for both avg (avg_gdd) and be (be_gdd) methods up until May 1st (doy: 121)
# calculate cumulative sum of CD (chill days) and CD^2 and CD^3 up until May 1st (doy: 121)
data <-
  historic_temperatures %>%
  mutate(avg_temp = ifelse(is.na(avg_temp), 0, avg_temp),
         gdd = ifelse(is.na(gdd), 0, gdd),
         year = parse_number(format(date, "%Y")),
         doy = date.to.DOY(format(date, "%Y/%m/%d"), format = 'yyyy/mm/dd')) %>%
  group_by(location, year) %>%
  nest() %>%
  left_join(cherry) %>%
  mutate(avg_gdd_sum = map(data, function(df) cumsum(df$avg_temp)[121]),
         avg_gdd2_sum = map(data, function(df) cumsum(df$avg_temp^2)[121]),
         avg_gdd3_sum = map(data, function(df) cumsum(df$avg_temp^3)[121]),
         be_gdd_sum = map(data, function(df) cumsum(df$gdd)[121]),
         be_gdd2_sum = map(data, function(df) cumsum(df$gdd^2)[121]),
         be_gdd3_sum = map(data, function(df) cumsum(df$gdd^3)[121]),
         cd_sum = map(data, function(df) cumsum(df$cd)[121])) %>%
  unnest(c(avg_gdd_sum, avg_gdd2_sum, avg_gdd3_sum, be_gdd_sum, be_gdd2_sum, be_gdd3_sum, cd_sum))
# View(data)

# find rows that contain NAs 
which(is.na(data), arr.ind=TRUE)
# row 130 and 55 contain NAs
# remove these rows from the dataset
data2 <- data[-c(55,130),]
# View(data2)

# remove columns 3-6 (not predictor variables)
data2 <- data2[,-c(1, 3:7)]
# View(data2)

# divide data into train/test sets
set.seed(130250)  

index <- createDataPartition(data2$bloom_doy, p=0.8, list=FALSE, times=1)
train_df <- data2[index,]
test_df <- data2[-index,]

y_train <- train_df$bloom_doy
x_train <- data.matrix(train_df[, c("avg_gdd_sum", "avg_gdd2_sum", "avg_gdd3_sum", "be_gdd_sum", "be_gdd2_sum", "be_gdd3_sum", "cd_sum")])
train <- cbind(x_train, y_train)

y_test <- test_df$bloom_doy
x_test <- data.matrix(test_df[, c("avg_gdd_sum", "avg_gdd2_sum", "avg_gdd3_sum", "be_gdd_sum", "be_gdd2_sum", "be_gdd3_sum", "cd_sum")])
test <- cbind(x_test, y_test)

# elastic net regression model fitting
control <- trainControl(method = "cv",
                        number = 5)

elastic_model <- train(y_train ~ .,
                       data = train,
                       method = "glmnet",
                       preProcess = c("center", "scale"),
                       tuneLength = 25,
                       trControl = control)

elastic_model

y_predicted <- predict(elastic_model, x_test)
y_predicted

# model summary
mod_summary <- function(actual, predicted) {
  sst <- sum((actual - mean(actual))^2)
  sse <- sum((predicted - actual)^2)
  rsq <- 1 - sse/sst
  rmse <- sqrt(mean((predicted - actual)^2))
  mae <- MAE(predicted, actual)
  output <- c(rsq, rmse, mae)
  print("R-squared & RMSE & MAE")
  return(output)
}

mod_summary(actual = y_test, predicted = y_predicted)

# model plots
plot(elastic_model, main = "Elastic Net Regression")
plot(varImp(elastic_model, scale = TRUE))

# find forecasted temperature predictions for the next decade from Accuweather.com
get_weather_table <- function(url)
  read_html(url) %>% 
  html_nodes("div.monthly-calendar") %>% 
  html_text2() %>%
  str_remove_all("°|Hist. Avg. ") %>%
  str_replace("N/A", "NA NA") %>%
  str_split(" ", simplify = TRUE) %>%
  parse_number() %>%
  matrix(ncol = 3, 
         byrow = TRUE,
         dimnames = list(NULL, c("day", "tmax", "tmin"))) %>%
  as_tibble() %>%
  filter(
    row_number() %in%
      (which(diff(day) < 0) %>% (function(x) if(length(x) == 1) seq(1, x[1], 1) else seq(x[1] + 1, x[2], 1))))

accu_temp <- function(year, city){
  tibble(
    base_url = paste("accuweather", year, "/", city ,"/", sep = ""),
    month = tolower(month.abb)[1:5],
    url = str_c(base_url, month, ".html")) %>%
    mutate(temp = map(url, get_weather_table)) %>%
    pull(temp) %>%
    reduce(bind_rows) %>%
    transmute(date = seq(as.Date(paste(year, "-01-01", sep = "")), as.Date(paste(year, "-05-31", sep = "")), 1),
              doy = date.to.DOY(format(date, "%Y/%m/%d"), format = 'yyyy/mm/dd'),
              year = parse_number(format(date, "%Y")),
              tmax = (tmax - 32) * 5/9,
              tmin = (tmin - 32) * 5/9,
              temp = (tmax + tmin) / 2)
}

setup <- function (doy, tmin, tmax) {
  avg_temp <- (tmax + tmin) / 2
  W <- (tmax - tmin)/2
  A <- asin((base_temp - avg_temp) / W)
  gdd <- case_when(
    ((tmax < base_temp) | (is.na(tmin) == TRUE) | (is.na(tmax) == TRUE)) ~ 0,
    (tmin >= base_temp) ~ avg_temp - base_temp,
    (tmin < base_temp) ~ ((W * cos(A)) - ((base_temp - avg_temp) * ((3.14/2)-A)))/3.14)
  cd <- case_when(
    (tmax < base_temp) ~ 1,
    .default = 0)
  variable_df <- as.data.frame(cbind(doy, avg_temp, gdd, cd))
  colnames(variable_df) <- c("doy", "avg_temp", "gdd", "cd")
  return(variable_df)
}

predictors <- function (df) {
  avg_gdd_sum = cumsum(df$avg_temp)[121]
  avg_gdd2_sum = cumsum(df$avg_temp^2)[121]
  avg_gdd3_sum = cumsum(df$avg_temp^3)[121]
  be_gdd_sum = cumsum(df$gdd)[121]
  be_gdd2_sum = cumsum(df$gdd^2)[121]
  be_gdd3_sum = cumsum(df$gdd^3)[121]
  cd_sum = cumsum(df$cd)[121]
  predictors_df <- as.data.frame(cbind(avg_gdd_sum, avg_gdd2_sum, avg_gdd3_sum, be_gdd_sum, be_gdd2_sum, be_gdd3_sum, cd_sum))
  colnames(predictors_df) <- c("avg_gdd_sum", "avg_gdd2_sum", "avg_gdd3_sum", "be_gdd_sum", "be_gdd2_sum", "be_gdd3_sum", "cd_sum")
  return(predictors_df)
}

# find temperature data & forecasted data for 2023
dc_2023_temp <- as.data.frame(accu_temp("2023", "washingtondc"))
kyoto_2023_temp <- as.data.frame(accu_temp("2023","kyoto"))
liestal_2023_temp <- as.data.frame(accu_temp("2023","liestal"))
vancouver_2023_temp <- as.data.frame(accu_temp("2023","vancouver"))

# impute missing day for kyoto data
kyoto_2023_temp[59, 4] <- (kyoto_2023_temp[58, 4] + kyoto_2023_temp[60, 4]) / 2
kyoto_2023_temp[59, 5] <- (kyoto_2023_temp[58, 5] + kyoto_2023_temp[60, 5]) / 2
kyoto_2023_temp[59, 6] <- (kyoto_2023_temp[59, 4] + kyoto_2023_temp[59, 5]) / 2

# find predictor values for 2023
dc_2023_pred <- round(predict(elastic_model, predictors(setup(dc_2023_temp$doy, dc_2023_temp$tmin, dc_2023_temp$tmax))))
kyoto_2023_pred <- round(predict(elastic_model, predictors(setup(kyoto_2023_temp$doy, kyoto_2023_temp$tmin, kyoto_2023_temp$tmax))))
liestal_2023_pred <- round(predict(elastic_model, predictors(setup(liestal_2023_temp$doy, liestal_2023_temp$tmin, liestal_2023_temp$tmax))))
vancouver_2023_pred <- round(predict(elastic_model, predictors(setup(vancouver_2023_temp$doy, vancouver_2023_temp$tmin, vancouver_2023_temp$tmax))))

pred_2023 <- cbind(kyoto_2023_pred, liestal_2023_pred, dc_2023_pred, vancouver_2023_pred)

# find historic temperature data for 2024-2032 (same historic temps used)
dc_historic_temp <- as.data.frame(accu_temp("2024", "washingtondc"))
kyoto_historic_temp <- as.data.frame(accu_temp("2024","kyoto"))
liestal_historic_temp <- as.data.frame(accu_temp("2024","liestal"))
vancouver_historic_temp <- as.data.frame(accu_temp("2024","vancouver"))

# find predictions for 2024-2032
dc_future_pred <- round(predict(elastic_model, predictors(setup(dc_historic_temp$doy, dc_historic_temp$tmin, dc_historic_temp$tmax))))
kyoto_future_pred <- round(predict(elastic_model, predictors(setup(kyoto_historic_temp$doy, kyoto_historic_temp$tmin, kyoto_historic_temp$tmax))))
liestal_future_pred <- round(predict(elastic_model, predictors(setup(liestal_historic_temp$doy, liestal_historic_temp$tmin, liestal_historic_temp$tmax))))
vancouver_future_pred <- round(predict(elastic_model, predictors(setup(vancouver_historic_temp$doy, vancouver_historic_temp$tmin, vancouver_historic_temp$tmax))))

pred_future <- cbind(kyoto_future_pred, liestal_future_pred, dc_future_pred, vancouver_future_pred)

# compile predictions
year <- c(2023:2032)
predictions <- rbind(pred_2023, pred_future)
predictions <- rbind(predictions, predictions[rep(2, 8),])
predictions <- as.data.frame(cbind(year, predictions))
rownames(predictions) <- NULL
colnames(predictions) <- c("year","kyoto","liestal","washingtondc","vancouver")
View(predictions)

# write to .csv file
write.csv(predictions, "C:\\Users\\madis\\Downloads\\predictions.csv", row.names=FALSE)

