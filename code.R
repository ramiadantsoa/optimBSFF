library(lubridate)
library(tidyverse)
library(patchwork)
library(lme4)
library(strucchange)

fao.yield.maize <- 1.81 # t per ha
fao.yield.soy <- 0.6 # t per ha

fao.price.maize <- 130.4 # USD per ha
fao.price.soy <- 352.4 # USD per ha

date.soil.collect <- dmy("27-04-2024")
date.soil.analyses <- dmy("10-07-2024")
date.soil.prep <- dmy("07-01-2025")
date.plot.create <- dmy("08-01-2025")
date.plot.frass0 <- dmy("10-01-2025")
date.plot.frass1 <- dmy("13-01-2025")
date.semi.maize <- dmy("14-01-2025")
date.semi.soy <- dmy("15-01-2025")
date.flood <- dmy("23-01-2025")
date.pyrifos1 <- dmy("02-02-2025")
date.pyrifos2 <- dmy("04-02-2025")
date.megalegion1 <- dmy("06-02-2025")
date.megalegion2 <- dmy("08-02-2025")

date.germination <- dmy("30-01-2025")
date.growth <- dmy("22-02-2025")
date.flower1 <- dmy("03-03-2025")
date.flower2 <- dmy("10-03-2025")
date.yield.soy1a <- dmy("15-04-2025")
date.yield.soy1b <- dmy("16-04-2025")
date.yield.soy2 <- dmy("17-04-2025")
date.yield.maize1 <- dmy("17-04-2025")
date.yield.maize2a <- dmy("19-05-2025")
date.yield.maize2b <- dmy("22-05-2025")

labsx <- labs(x = expression("Fertilizer application rate (t ha"^-1 * ")"))
labsy_yield <- labs(y = expression("Total yield (t ha"^-1 * ")"))
col.dol.1 <- "blue"
col.dol.0 <- "red"
shape.flood.Y <- 19
shape.flood.N <- 1
my_theme <- theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
)

dma <- read.csv("data/maize_plot_level.csv")
dm <- read.csv("data/maize_plant_level.csv")
dsa <- read.csv("data/soy_plot_level.csv")
ds <- read.csv("data/soy_plant_level.csv")


dma$Dolomite <- as.character(dma$Dolomite)
dm$Dolomite <- as.character(dm$Dolomite)
dsa$Dolomite <- as.character(dsa$Dolomite)
ds$Dolomite <- as.character(ds$Dolomite)
# number of individuals for maize and soybean
nm <- 30
ns <- 39

# outlier plots for maize and soybean
out.plot.m <- c(39, 38, 64, 59)
out.plot.s <- c(94, 93)

# flooded plot from all plots
# flooded.maize <- dma |> filter(flooded == "Y") |> count()
# flooded.soy <- dsa |> filter(flooded == "Y") |> count()

# unflooded maize and soy
unflooded.maize <- dma$plotID |> unique() |> length()
unflooded.soy <- dsa$plotID |> unique() |> length()

# prediction using logit model for glmer
f.mod <- function(x, mod, n) {
    aa <- fixef(mod)["Doses"]
    bb <- fixef(mod)["(Intercept)"]
    y = n * exp( aa * x + bb) / (1 + exp( aa * x + bb))
    y
}

# prediction using logit model for glm
f.mod2 <- function(x, mod, n) {
    aa <- coef(mod)["Doses"]
    bb <- coef(mod)["(Intercept)"]
    y = n * exp( aa * x + bb) / (1 + exp( aa * x + bb))
    y
}

# prediction using lm model for glm
f.mod.lm <- function(x, mod) {
    aa <- coef(mod)["Doses"]
    bb <- coef(mod)["(Intercept)"]
    y =  aa * x + bb
    y
}

# extraction slope from linear model
lm_slope <- function(x, y){
    mod <- lm(y~x)
    slope <- coef(mod)[2]
    slope
}

# prediction based on asymptotic model
f.Asym <- function(x, mod) {
    cc <- coef(mod)
    asym <- cc["Asym"]
    r0 <- cc["R0"]
    lrc <- cc["lrc"]
    y  <- asym + (r0 - asym) * exp(- exp(lrc) * x)
    y
}

print.date <- function(x){
    format(x, "%e %B, %Y")
}

print.pvalue <- function(value){
    if(value <= 0.001){
        res <- " < 0.001"
    }else {
        res <- paste0(" = ", round(value,3))
    }
    return(res)
}


#####################
#### Germination ####
#####################

d.ger.maize <- dma |>
    select(Doses, nGer0, plotID) |>
    mutate(
        germination = cbind(nGer0, nonGer = nm - nGer0))

mod.ger.maize <- glmer(germination ~ Doses + (1|plotID), family = binomial, data = d.ger.maize)

d.ger.soy <- dsa |>
    select(Doses,  nGer0, plotID) |>
    mutate(
        germination = cbind(nGer0, nonGer = ns - nGer0))

mod.ger.soy <- glmer(germination ~ Doses + (1|plotID), family = binomial, data = d.ger.soy)


f.ger.maize <- dma |>
    filter(!plotID %in% out.plot.m) |>
    ggplot( aes(x = Doses, y = nGer0 / nm)) +
        geom_jitter(height = 0, shape = 1) +
        geom_function(fun = ~f.mod(.x, mod.ger.maize, 1), lwd = 1.5) +
        ylim(c(0, 1)) +
        theme_bw() +
        labsx +
        labs(y = "Germination rate") +
        my_theme

f.ger.soy <- ggplot(dsa, aes(x = Doses, y = nGer0 / ns)) +
        geom_jitter(height = 0, shape = 1) +
            geom_function(fun = ~f.mod(.x, mod.ger.soy, 1), lwd = 1.5) +
            ylim(c(0, 1)) +
            theme_bw() +
            labsx +
            labs(y = "Germination rate") +
            my_theme

fig.ger <- f.ger.maize + f.ger.soy +
    plot_annotation(tag_levels = 'A')


#####################
#### Yield individual  ####
#####################

d.yieldI.maize <- dm |>
    filter(!is.na(weightAll), !plotID %in% out.plot.m) |>
    group_by(Doses, plotID) |>
    summarise(m.wAll = mean(weightAll))

mod.yieldI.maize <- nls(m.wAll ~ SSasymp(Doses, Asym, R0, lrc), data = d.yieldI.maize)

d.yieldI.soy <- ds |>
    filter(!is.na(weightAll) ) |>
    group_by(Doses, plotID) |>
    summarise(m.wAll = mean(weightAll))

mod.yieldI.soy <- nls(m.wAll ~ SSasymp(Doses, Asym, R0, lrc), data = d.yieldI.soy)

d.slope.maize <- dm |>
    filter(!is.na(heightLeaf) & !is.na(weightAll)) |>
    group_by(Doses, plotID) |>
    summarise(slope = lm_slope(heightLeaf, weightAll))

mod.slope.maize <- lm(slope ~ Doses, data = d.slope.maize)

d.slope.soy <- ds |>
    filter(!is.na(cover) & !is.na(weightAll)) |>
    group_by(Doses, plotID) |>
    summarise(slope = lm_slope(cover, weightAll))

mod.slope.soy <- lm(slope ~ Doses , data = d.slope.soy)



f.yieldI.maize <- d.yieldI.maize |>
    ggplot(aes(x = Doses, y = m.wAll)) +
        geom_point(shape = 1) +
        geom_function(fun = ~f.Asym(.x, mod.yieldI.maize), lwd = 1.5) +
        theme_bw() +
        labsx +
        labs(y = "Mean yield per individual (g)") +
        my_theme

f.yieldI.soy <- d.yieldI.soy |>
    ggplot(aes(x = Doses, y = m.wAll)) +
        geom_point(shape = 1) +
        geom_function(fun = ~f.Asym(.x, mod.yieldI.soy), lwd = 1.5) +
        theme_bw() +
        labsx +
        labs(y = "Mean yield per individual (g)") +
        my_theme

f.slope.maize <- d.slope.maize |>
    ggplot(aes(x = Doses, y = slope)) +
        geom_jitter(shape = 1) +
        geom_smooth(method = "lm", se = FALSE, col = "black", lwd = 1.5) +
        theme_bw() +
        labsx +
        labs(y = "Slope of lm(yield ~ height)") +
        my_theme

f.slope.soy <- d.slope.soy  |>
    ggplot(aes(x = Doses, y = slope)) +
        geom_jitter(shape = 1) +
        geom_smooth(method = "lm", se = FALSE, col = "black", lwd = 1.5) +
        theme_bw() +
        labsx +
        labs(y = "Slope of lm(yield ~ cover)") +
        my_theme

fig.yield.individual <- (f.yieldI.maize + f.yieldI.soy) / (f.slope.maize + f.slope.soy) +
    plot_annotation(tag_levels = 'A')


#####################
#### Yield Total  ####
#####################
d.yieldA.maize <-  dma |>
    filter(!plotID %in% out.plot.m)
mod.yieldA.maize <- nls(yield ~ SSasymp(Doses, Asym, R0, lrc ), data = d.yieldA.maize)

d.yieldA.soy <-  dsa |>
    filter(!plotID %in% out.plot.s)
mod.yieldA.soy <- nls(yield ~ SSasymp(Doses, Asym, R0, lrc ), data = d.yieldA.soy)

Y0_maize <- f.Asym(0, mod.yieldA.maize)
Y0_soy <- f.Asym(0, mod.yieldA.soy)

Ylocal_maize <- f.Asym(1.2, mod.yieldA.maize)
Ylocal_soy <- f.Asym(1.2, mod.yieldA.soy)

Yfofifa_maize <- f.Asym(3.6, mod.yieldA.maize)
Yfofifa_soy <- f.Asym(3.6, mod.yieldA.soy)


# AE_maize1 <- (f.Asym(0.04, mod.yieldA.maize) - Y0_maize) / (0.04 * 0.0346)
AE_maize1 <- (f.Asym(0.4, mod.yieldA.maize) - Y0_maize) / (0.4 * 0.0346)
AE_maize2 <- (f.Asym(30, mod.yieldA.maize) - Y0_maize) / (30 * 0.0346)

relAE_maize <- (AE_maize1 - AE_maize2) / AE_maize1 * 100

AE_soy1 <- (f.Asym(0.4, mod.yieldA.soy) - Y0_soy) / (0.4 * 0.0346)
AE_soy2 <- (f.Asym(30, mod.yieldA.soy) - Y0_soy) / (30 * 0.0346)

relAE_soy <- (AE_soy1 - AE_soy2) / AE_soy1 * 100

# f_maize <- function(x){
#     y <- (f.Asym(x, mod.yieldA.maize) * (f.Asym(x, mod.yieldA.maize) - Y0_maize) / (x * 0.0346))
#     y
# }

# optim_maize <- optimize(f = f_maize, interval = c(0, 30), maximum = TRUE)$maximum
# optim_yield_maize <- f.Asym(optim_maize, mod.yieldA.maize)

# f_soy <- function(x){
#     y <- (f.Asym(x, mod.yieldA.soy) * (f.Asym(x, mod.yieldA.soy) - Y0_soy) / (x * 0.0346))
#     y
# }

# optim_soy <- optimize(f = f_soy, interval = c(0, 30), maximum = TRUE)$maximum
# optim_yield_soy <- f.Asym(optim_soy, mod.yieldA.soy)

########### temp remove this after #########
optim_maize <- 123
optim_soy <- 234

f.yieldA.maize <- ggplot(d.yieldA.maize, aes(x = Doses, y = yield)) +
    geom_point(shape = 1) +
    geom_function(fun = ~f.Asym(.x, mod.yieldA.maize), lwd = 1.5) +
    theme_bw() +
    labsx +
    labs(y = "Total yield (T/ha)") +
    my_theme

f.yieldA.soy <- ggplot(d.yieldA.soy, aes(x = Doses, y = yield)) +
    geom_point(shape = 1) +
    geom_function(fun = ~f.Asym(.x, mod.yieldA.soy), lwd = 1.5) +
    theme_bw() +
    labsx +
    labs(y = "Total yield (T/ha)") +
    my_theme

f.agroE <- ggplot(d.yieldA.maize, aes(x = Doses, y = yield)) +
    geom_function(
        fun = ~((f.Asym(.x, mod.yieldA.maize) - Y0_maize) / (.x * 0.0346)),
        lwd = 1.5,
        aes(linetype = "Maize")) +
    geom_function(
        fun = ~((f.Asym(.x, mod.yieldA.soy) - Y0_soy) / (.x * 0.0346)),
        lwd = 1.5,
        aes(linetype = "Soy")) +
    theme_bw() +
    labsx +
    labs(y = "NAE (kg / kg N)", linetype = "Crop") +
    scale_linetype_manual(values = c(Maize = 1, Soy = 2)) +
    theme(
        legend.position = c(0.8, 0.85),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.key.width = unit(3, "lines")
    ) +
    my_theme

# f.optim <- ggplot(d.yieldA.maize, aes(x = Doses, y = yield)) +
#     geom_function(
#         fun = ~(f.Asym(.x, mod.yieldA.maize) * (f.Asym(.x, mod.yieldA.maize) - Y0_maize) / (.x * 0.0346)) / 80,
#         lwd = 1.5,
#         aes(linetype = "Maize")) +
#     geom_function(
#         fun = ~(f.Asym(.x, mod.yieldA.soy) * (f.Asym(.x, mod.yieldA.soy) - Y0_soy) / (.x * 0.0346)) / 10,
#         lwd = 1.5,
#         aes(linetype = "Soy")) +
#     theme_bw() +
#     labsx +
#     labs(y = "NAE-adjusted yield", linetype = "Crop") +
#     scale_linetype_manual(values = c(Maize = 1, Soy = 2)) +
#     theme(
#         legend.position = c(0.8, 0.85),
#         legend.background = element_rect(fill = "white", color = "black"),
#         legend.key.width = unit(3, "lines")
#     ) +
#     my_theme

fig.yield.total <- f.yieldA.maize + f.yieldA.soy + f.agroE +
    plot_annotation(tag_levels = 'A') +
    plot_layout(axis_titles = "collect")


#####################
#### MVCR from Mathematica ####
#####################

# using 1.2 t/ha of frass
frassP.maize.local.1 <- 52.59
frassP.maize.local.2 <- 35.06
frassP.maize.local.3 <- 26.30
frassP.maize.fofifa.1 <- 46.20
frassP.maize.fofifa.2 <- 30.80
frassP.maize.fofifa.3 <- 23.10

mvcr.maize.price.soy.local.1 <- 2.32
mvcr.maize.price.soy.fofifa.1 <- 1.92

mvcr.maize.price.soy.local.3 <- 5.65
mvcr.maize.price.soy.fofifa.3 <- 4.85

# using 3 * 1.2 t/ha of frass
frassP.soy.local.1 <- 31.61
frassP.soy.local.2 <- 21.07
frassP.soy.local.3 <- 15.80
frassP.soy.fofifa.1 <- 27.75
frassP.soy.fofifa.2 <- 18.50
frassP.soy.fofifa.3 <- 13.87
