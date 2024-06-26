# txt_mod_fixed <- function(model_no = 1L) {
#
#   if (model_no == 1) {
#     mod <- "
#       eta1 =~ 0.80*y1 + 0.70*y2 + 0.47*y3 + 0.38*y4 + 0.34*y5
#
#       y1 | -1.43 * t1
#       y2 | -0.55 * t1
#       y3 | -0.13 * t1
#       y4 | -0.72 * t1
#       y5 | -1.13 * t1
#     "
#   }
#   if (model_no == 2) {
#     mod <- "
#       eta1 =~ 0.80*y1 + 0.70*y2 + 0.47*y3 + 0.38*y4 + 0.34*y5 +
#               0.80*y6 + 0.70*y7 + 0.47*y8
#
#       y1 | -1.43 * t1
#       y2 | -0.55 * t1
#       y3 | -0.13 * t1
#       y4 | -0.72 * t1
#       y5 | -1.13 * t1
#       y6 | -1.43 * t1
#       y7 | -0.55 * t1
#       y8 | -0.13 * t1
#     "
#   }
#   if (model_no == 3) {
#     mod <- "
#       eta1 =~ 0.80*y1  + 0.70*y2  + 0.47*y3  + 0.38*y4  + 0.34*y5 +
#               0.80*y6  + 0.70*y7  + 0.47*y8  + 0.38*y9  + 0.34*y10 +
#               0.80*y11 + 0.70*y12 + 0.47*y13 + 0.38*y14 + 0.34*y15 +
#
#       y1  | -1.43 * t1
#       y2  | -0.55 * t1
#       y3  | -0.13 * t1
#       y4  | -0.72 * t1
#       y5  | -1.13 * t1
#       y6  | -1.43 * t1
#       y7  | -0.55 * t1
#       y8  | -0.13 * t1
#       y9  | -0.72 * t1
#       y10 | -1.13 * t1
#       y11 | -1.43 * t1
#       y12 | -0.55 * t1
#       y13 | -0.13 * t1
#       y14 | -0.72 * t1
#       y15 | -1.13 * t1
#     "
#   }
#   if (model_no == 4) {
#     mod <- "
#       eta1 =~ 0.80*y1  + 0.70*y2  + 0.47*y3  + 0.38*y4  + 0.34*y5 +
#       eta2 =~ 0.80*y6  + 0.70*y7  + 0.47*y8  + 0.38*y9  + 0.34*y10
#
#       y1  | -1.43 * t1
#       y2  | -0.55 * t1
#       y3  | -0.13 * t1
#       y4  | -0.72 * t1
#       y5  | -1.13 * t1
#       y6  | -1.43 * t1
#       y7  | -0.55 * t1
#       y8  | -0.13 * t1
#       y9  | -0.72 * t1
#       y10 | -1.13 * t1
#
#       eta1 ~~ 0.30*eta2
#     "
#   }
#   if (model_no == 5) {
#     mod <- "
#       eta1 =~ 0.80*y1  + 0.70*y2  + 0.47*y3  + 0.38*y4  + 0.34*y5
#       eta2 =~ 0.80*y6  + 0.70*y7  + 0.47*y8  + 0.38*y9  + 0.34*y10
#       eta3 =~ 0.80*y11 + 0.70*y12 + 0.47*y13 + 0.38*y14 + 0.34*y15
#
#       y1  | -1.43 * t1
#       y2  | -0.55 * t1
#       y3  | -0.13 * t1
#       y4  | -0.72 * t1
#       y5  | -1.13 * t1
#       y6  | -1.43 * t1
#       y7  | -0.55 * t1
#       y8  | -0.13 * t1
#       y9  | -0.72 * t1
#       y10 | -1.13 * t1
#       y11 | -1.43 * t1
#       y12 | -0.55 * t1
#       y13 | -0.13 * t1
#       y14 | -0.72 * t1
#       y15 | -1.13 * t1
#
#       eta1 ~~ 0.2*eta2
#       eta1 ~~ 0.3*eta3
#       eta2 ~~ 0.4*eta3
#     "
#   }
#
#   mod
# }
