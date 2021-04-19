# WGIF-and-GIF
This project implements Guided image filter, Weighted guided image filter, and SD filter in Julia with its extensions.
## Usage of Filters
src/gif.jl contains functionn gif(...) for guided image filter.

src/wgif.jl contains function wgif(...) and its extensions wgif_1(...), wgif_2(...), wgif_3(...), wgif_4(...) for weighted guided image filter.

src/rgif.jl contains function rgif(...) for rolling guided image filter.

src/sdfilter.jl contains function sdfilter(...), asd(...), agsd(...) for SD filter, adaptive SD filter, and augmented guidence SD filter. 

src/main.jl contains the default parameters of each filters and its application.
