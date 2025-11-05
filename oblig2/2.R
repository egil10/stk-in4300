library(tidyverse)
library(mlbench)
library(class)
library(mgcv)
library(rpart)
library(rpart.plot)
library(ipred)
library(randomForest)
set.seed(4300)
theme_set(theme_bw())

# ------------------------------------------------------------
# Problem 2 — Classification on PimaIndiansDiabetes
# ------------------------------------------------------------

data("PimaIndiansDiabetes")
df <- as_tibble(PimaIndiansDiabetes)

# ---------------------------
# (a) Train/test split + kNN
# ---------------------------

# Stratified split: keep class balance
ix_pos <- which(df$diabetes == "pos")
ix_neg <- which(df$diabetes == "neg")

tr_ix <- c(
  sample(ix_pos, round(2 * length(ix_pos) / 3)),
  sample(ix_neg, round(2 * length(ix_neg) / 3))
)

tr <- df[tr_ix, ]
te <- df[-tr_ix, ]

# Standardize predictors using training mean/sd
Xtr <- scale(dplyr::select(tr, -diabetes))
ctr <- attr(Xtr, "scaled:center")
str <- attr(Xtr, "scaled:scale")

Xte <- scale(dplyr::select(te, -diabetes), center = ctr, scale = str)

ytr <- tr$diabetes
yte <- te$diabetes

kvals <- 1:30
cv5   <- rep(NA, length(kvals))
loocv <- rep(NA, length(kvals))
te_err <- rep(NA, length(kvals))

# Stratified 5-fold CV
folds <- 5
pos_f <- split(sample(which(ytr == "pos")), rep(1:folds, length.out = sum(ytr == "pos")))
neg_f <- split(sample(which(ytr == "neg")), rep(1:folds, length.out = sum(ytr == "neg")))

for(i in seq_along(kvals)){
  k <- kvals[i]
  
  # 5-fold CV error
  cv5[i] <- mean(sapply(1:folds, function(f){
    tr_id <- c(unlist(pos_f[-f]), unlist(neg_f[-f]))
    va_id <- c(unlist(pos_f[f]),  unlist(neg_f[f]))
    pr <- knn(train = Xtr[tr_id,], test = Xtr[va_id,], cl = ytr[tr_id], k = k)
    mean(pr != ytr[va_id])
  }))
  
  # LOOCV error
  loocv[i] <- mean(knn.cv(Xtr, ytr, k = k) != ytr)
  
  # Test error
  te_err[i] <- mean(knn(train = Xtr, test = Xte, cl = ytr, k = k) != yte)
}

# Plot all three error curves
p <- tibble(k = kvals, cv5 = cv5, loocv = loocv, test = te_err) %>%
  pivot_longer(-k, names_to = "type", values_to = "err") %>%
  ggplot(aes(k, err, linetype = type)) + geom_line()

ggsave("plots/plot2a.pdf", p, width = 10, height = 6)

cat("2a_cv5_min", min(cv5),
    "2a_loocv_min", min(loocv),
    "2a_test_at_k_cv5", te_err[which.min(cv5)], "\n")


# ---------------------------
# (b) GAM model
# ---------------------------

g_full <- gam(
  diabetes ~ s(pregnant) + s(glucose) + s(pressure) + s(triceps) +
    s(insulin) + s(mass) + s(pedigree) + s(age),
  family = binomial,
  data = tr,
  method = "REML",
  select = TRUE
)

# Classify based on predicted probability > 0.5
pr_b <- ifelse(predict(g_full, te, type = "response") > 0.5, "pos", "neg")
te_gam <- mean(pr_b != yte)

cat("2b_gam_test", te_gam, "\n")


# ---------------------------
# (c) Trees, Bagging, Random Forest
# ---------------------------

# CART tree + prune using optimal CP
tree <- rpart(diabetes ~ ., data = tr, method = "class", cp = 0.001)
cp_opt <- tree$cptable[which.min(tree$cptable[, "xerror"]), "CP"]
tree <- prune(tree, cp = cp_opt)

p_tr <- predict(tree, tr, type = "class")
p_te <- predict(tree, te, type = "class")

tr_tree <- mean(p_tr != ytr)
te_tree <- mean(p_te != yte)

# Save tree plot
pdf("plots/plot2c_tree.pdf", width = 10, height = 6)
rpart.plot(tree)
dev.off()

# Bagging
bag <- bagging(diabetes ~ ., data = tr, nbagg = 200)
tr_bag <- mean(predict(bag, tr, type = "class") != ytr)
te_bag <- mean(predict(bag, te, type = "class") != yte)

# Random Forest
rf <- randomForest(diabetes ~ ., data = tr, ntree = 500)
tr_rf <- mean(predict(rf, tr) != ytr)
te_rf <- mean(predict(rf, te) != yte)

cat("2c_tree", tr_tree, te_tree,
    "2c_bag", tr_bag, te_bag,
    "2c_rf", tr_rf, te_rf, "\n")


# ---------------------------
# (d) Best model by test error
# ---------------------------

cat("2d_best_by_test",
    c("tree" = te_tree,
      "bag" = te_bag,
      "rf"  = te_rf,
      "gam" = te_gam)[ which.min(c(te_tree, te_bag, te_rf, te_gam)) ],
    "\n")


# ---------------------------
# (e) Repeat analysis on cleaned dataset
# ---------------------------

data("PimaIndiansDiabetes2")
df2 <- as_tibble(PimaIndiansDiabetes2) %>% drop_na()

# Stratified split again
ix_pos2 <- which(df2$diabetes == "pos")
ix_neg2 <- which(df2$diabetes == "neg")

tr_ix2 <- c(
  sample(ix_pos2, round(2 * length(ix_pos2) / 3)),
  sample(ix_neg2, round(2 * length(ix_neg2) / 3))
)

tr2 <- df2[tr_ix2, ]
te2 <- df2[-tr_ix2, ]

# Standardize
Xtr2 <- scale(dplyr::select(tr2, -diabetes))
ctr2 <- attr(Xtr2, "scaled:center")
str2 <- attr(Xtr2, "scaled:scale")

Xte2 <- scale(dplyr::select(te2, -diabetes), center = ctr2, scale = str2)

ytr2 <- tr2$diabetes
yte2 <- te2$diabetes

# Use best k from (a)
k <- kvals[which.min(cv5)]

te_knn2 <- mean(knn(train = Xtr2, test = Xte2, cl = ytr2, k = k) != yte2)

# GAM on cleaned data
g2 <- gam(
  diabetes ~ s(pregnant) + s(glucose) + s(pressure) + s(triceps) +
    s(insulin) + s(mass) + s(pedigree) + s(age),
  family = binomial,
  data = tr2,
  method = "REML",
  select = TRUE
)

te_gam2 <- mean(ifelse(predict(g2, te2, type = "response") > 0.5, "pos", "neg") != yte2)

# Tree
tree2 <- rpart(diabetes ~ ., data = tr2, method = "class", cp = 0.001)
cp_opt2 <- tree2$cptable[which.min(tree2$cptable[, "xerror"]), "CP"]
tree2 <- prune(tree2, cp = cp_opt2)
te_tree2 <- mean(predict(tree2, te2, type = "class") != yte2)

# Bagging
bag2 <- bagging(diabetes ~ ., data = tr2, nbagg = 200)
te_bag2 <- mean(predict(bag2, te2, type = "class") != yte2)

# Random Forest
rf2 <- randomForest(diabetes ~ ., data = tr2, ntree = 500)
te_rf2 <- mean(predict(rf2, te2) != yte2)

cat("2e_knn", te_knn2,
    "2e_gam", te_gam2,
    "2e_tree", te_tree2,
    "2e_bag", te_bag2,
    "2e_rf", te_rf2, "\n")
