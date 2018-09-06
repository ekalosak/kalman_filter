# Author: Eric Kalosa-Kenyon
# Analyzing poll data with a Kalman filter

## @knitr setup
library(ggplot2)
library(scales)

# https://raw.githubusercontent.com/fivethirtyeight/data/master/pollster-ratings/2016/raw-polls.csv
d = read.csv("raw538polls.csv", header=T)
n = dim(d)[1]
k = dim(d)[2]

# Data consists of 7977 observations going back to 1998 up to 2016
# There are 24 columns - id, race, loc, pollster, estimates, biases, and errors.

# TODO:
# Dropping the years with fewer than 50 observations for now
df_plot_year = as.data.frame(table(d$year))

## @knitr plot_year
df_plot_year = as.data.frame(table(d$year))
colnames(df_plot_year) = c("Year", "CountObsv")
p1 = ggplot(
            df_plot_year,
            aes(
                x=factor(1),
                y=CountObsv,
                fill=Year
                )
            ) +
    geom_bar(width=1, stat="identity") +
    coord_polar(theta="y", start=0) +
    geom_text(
                aes(
                    y = CountObsv/2 + c(0, cumsum(CountObsv)[-length(CountObsv)]),
                    # y = c(0, cumsum(CountObsv)[-length(CountObsv)]),
                    label = paste(percent(CountObsv/n), Year)
                    ),
            size=4,
            # nudge_x=1/2,
            # position=position_jitter(width=0.02)
            ) +
    labs(title="Observations per year", x="", y="")
