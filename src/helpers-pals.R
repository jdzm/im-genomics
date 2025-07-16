# My palettes

## Heatmap 
type_colors <- c("milk"   = "#045a8d",
                 "brine"   = "#2b8cbe",
                 "cheese" = "#74a9cf",
                 "crust" = "#bdc9e1",
                 "wall" = "#f1eef6")

test_colors <- c("control" = "#DFC27D", 
                 "none" = "#8AAC86",
                 "test"   = "#35978F")

time_colors <- c("early" = "#FEE0D2", "mid" = "#FC9272",
                 "late"   = "#DE2D26")

cluster_colors = c("1"="#F7F7F7",
                   "2"="#D9D9D9",
                   "3"="#BDBDBD",
                   "4"="#969696",
                   "5"="#636363",
                   "6"="#252525")

map_color = colorRampPalette(c("white","#238b45"))(100)

### CIPAIS
qual_colors <- c("Buono"    = "#003f5c", 
                 "Discreto" = "#7a5195",
                 "Mediocre" = "#ef5675",
                 "Cattivo"  = "#ffa600")

site_colors <- c("B1"   = "#003f5c","B2"= "#277191","B3"= "#4aa8c8","B4"= "#6fe2ff",
                 "D1"   = "#7a5195","D2"= "#a375b7","D3"= "#cc9cda","D4"= "#f7c4ff",
                 "M1"   = "#ef5675","M2"= "#f67890","M3"= "#fc97ab","M4"= "#ffb5c4",
                 "C1"   = "#ffa600","C2"= "#ffb33e","C3"= "#ffc063","C4"= "#ffcc84")



room_descs = c("Atelier architettura",
               "Entrata, inizio rampa",
               "Ufficio aperto ISAAC",
               "Ufficio aperto IMC",
               "Lato Stazione",
               "Mensa", 
               "Ufficio aperto IM", 
               "Aula 60 persone",
               "Aula-Atelier architettura")

# room_descs = setNames(palette.colors (n=9), room_descs)
room_descs = setNames(hcl.colors (n=9, palette = "Spectral"), room_descs)




