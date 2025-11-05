library(tidyverse)
library(MASS)
library(rpart)
library(rpart.plot)
library(glmnet)
library(mgcv)
set.seed(4300)
theme_set(theme_bw())

# ========================================================
# Problem 1: QSAR Aquatic Toxicity Regression
# ========================================================

# --------------------------------------------------------
# Load and prepare data
# --------------------------------------------------------
df <- read.csv("data/qsar_aquatic_toxicity.csv", sep=";", header=FALSE) %>%
  as_tibble()

colnames(df) <- c("TPSA","SAacc","H050","MLOGP","RDCHI","GATS1p","nN","C040","LC50")

# Train/test split (≈ 2/3 training, 1/3 test)
i  <- sample(1:nrow(df), round(2 * nrow(df) / 3))
tr <- df[i, ]
te <- df[-i, ]

# --------------------------------------------------------
# (a) Linear models: raw counts vs. dummy-encoded counts
# --------------------------------------------------------

# Model with raw counts included as linear terms
m1 <- lm(LC50 ~ ., data = tr)

# Dummy model: use 0/1 indicators for count variables
m2 <- lm(LC50 ~ TPSA + SAacc + MLOGP + RDCHI + GATS1p +
           I(H050 > 0) + I(nN > 0) + I(C040 > 0),
         data = tr)

# Predictions and train/test MSE
p1 <- predict(m1, tr);  q1 <- predict(m1, te)
p2 <- predict(m2, tr);  q2 <- predict(m2, te)

cat("train_mse_m1", mean((tr$LC50 - p1)^2),
    "test_mse_m1",  mean((te$LC50 - q1)^2), "\n")
cat("train_mse_m2", mean((tr$LC50 - p2)^2),
    "test_mse_m2",  mean((te$LC50 - q2)^2), "\n")

# Scatterplot of observed vs predicted
p <- ggplot(tibble(obs = te$LC50, pred = q1),
            aes(obs, pred)) +
  geom_point() +
  geom_abline()
ggsave("plots/plot1a.pdf", p, width=10, height=6)

# --------------------------------------------------------
# (b) Repeat experiment 200 times to estimate test error distribution
# --------------------------------------------------------
e1 <- numeric(200)  # raw counts model
e2 <- numeric(200)  # dummy model

for (k in 1:200) {
  i  <- sample(1:nrow(df), round(2 * nrow(df) / 3))
  tr <- df[i, ]; te <- df[-i, ]
  
  m1 <- lm(LC50 ~ ., data=tr)
  m2 <- lm(LC50 ~ TPSA + SAacc + MLOGP + RDCHI + GATS1p +
             I(H050 > 0) + I(nN > 0) + I(C040 > 0),
           data=tr)
  
  e1[k] <- mean((te$LC50 - predict(m1, te))^2)
  e2[k] <- mean((te$LC50 - predict(m2, te))^2)
}

# Boxplot of empirical MSE distributions
t <- tibble(model = rep(c("raw","dummy"), each=200),
            mse   = c(e1, e2))

p <- ggplot(t, aes(model, mse)) + geom_boxplot()
ggsave("plots/plot1b.pdf", p, width=10, height=6)

cat("mean_test_mse_raw", mean(e1), "\n")
cat("mean_test_mse_dummy", mean(e2), "\n")

# --------------------------------------------------------
# (c) Variable selection: forward/backward with AIC and BIC
# --------------------------------------------------------
full <- lm(LC50 ~ ., data=tr)
null <- lm(LC50 ~ 1, data=tr)

b_aic <- step(full, direction="backward", k=2, trace=0)
f_aic <- step(null, scope=formula(full), direction="forward", k=2, trace=0)
b_bic <- step(full, direction="backward", k=log(nrow(tr)), trace=0)
f_bic <- step(null, scope=formula(full), direction="forward", k=log(nrow(tr)), trace=0)

cat("backward_AIC:",  deparse(formula(b_aic)), "\n")
cat("forward_AIC :",  deparse(formula(f_aic)), "\n")
cat("backward_BIC:", deparse(formula(b_bic)), "\n")
cat("forward_BIC :", deparse(formula(f_bic)), "\n")

# --------------------------------------------------------
# (d) Ridge regression: CV vs bootstrap OOB estimate
# --------------------------------------------------------
xtr <- as.matrix(dplyr::select(tr, -LC50))
ytr <- tr$LC50
xte <- as.matrix(dplyr::select(te, -LC50))
yte <- te$LC50

lambda_grid <- 10^seq(-4, 4, length=100)

# Cross-validation
cv <- cv.glmnet(xtr, ytr, alpha=0, lambda=lambda_grid, standardize=TRUE)
lam_cv <- cv$lambda.min

tr_cv <- mean((ytr - predict(cv$glmnet.fit, s=lam_cv, newx=xtr))^2)
te_cv <- mean((yte - predict(cv$glmnet.fit, s=lam_cv, newx=xte))^2)

# Bootstrap OOB
B <- 200
oob_mse <- matrix(NA, B, length(lambda_grid))

for(b in 1:B){
  idx <- sample(seq_len(nrow(tr)), replace=TRUE)
  oob <- setdiff(seq_len(nrow(tr)), unique(idx))
  
  fit <- glmnet(xtr[idx,], ytr[idx], alpha=0, lambda=lambda_grid, standardize=TRUE)
  
  if(length(oob) > 5){
    pred_oob <- predict(fit, newx=xtr[oob,])
    oob_mse[b,] <- colMeans((ytr[oob] - pred_oob)^2)
  }
}

boot_curve <- colMeans(oob_mse, na.rm=TRUE)
lam_boot <- lambda_grid[which.min(boot_curve)]

# Save CV vs bootstrap plot
p <- tibble(log_lambda=log(lambda_grid),
            cv=cv$cvm[match(lambda_grid, cv$lambda)],
            boot=boot_curve) %>%
  pivot_longer(-log_lambda, names_to="method", values_to="mse") %>%
  ggplot(aes(log_lambda, mse, linetype=method)) +
  geom_line()

ggsave("plots/plot1d.pdf", p, width=10, height=6)

# --------------------------------------------------------
# (e) GAM models with different smoothness
# --------------------------------------------------------
g1 <- gam(LC50 ~ s(TPSA,k=4) + s(SAacc,k=4) + s(H050,k=3) +
            s(MLOGP,k=4) + s(RDCHI,k=4) + s(GATS1p,k=4) +
            s(nN,k=3) + s(C040,k=3),
          data=tr, method="REML")

g2 <- gam(LC50 ~ s(TPSA,k=7) + s(SAacc,k=7) + s(H050,k=3) +
            s(MLOGP,k=7) + s(RDCHI,k=7) + s(GATS1p,k=7) +
            s(nN,k=3) + s(C040,k=3),
          data=tr, method="REML")

tr_g1 <- mean((tr$LC50 - predict(g1,tr))^2)
te_g1 <- mean((te$LC50 - predict(g1,te))^2)
tr_g2 <- mean((tr$LC50 - predict(g2,tr))^2)
te_g2 <- mean((te$LC50 - predict(g2,te))^2)

# --------------------------------------------------------
# (f) Regression tree + cost-complexity pruning
# --------------------------------------------------------
tree0 <- rpart(LC50 ~ ., data=tr, method="anova", cp=0.001)
cp_opt <- tree0$cptable[which.min(tree0$cptable[,"xerror"]), "CP"]

tree <- prune(tree0, cp=cp_opt)

tr_tree <- mean((tr$LC50 - predict(tree,tr))^2)
te_tree <- mean((te$LC50 - predict(tree,te))^2)

pdf("plots/plot1f_tree.pdf", width=10, height=6); rpart.plot(tree); dev.off()
pdf("plots/plot1f_cp.pdf",   width=10, height=6); plotcp(tree0); dev.off()

# --------------------------------------------------------
# (g) Compare all models
# --------------------------------------------------------
res <- tibble(
  model=c("lm_raw","lm_dummy","step_back_AIC","step_fwd_AIC",
          "step_back_BIC","step_fwd_BIC","ridge_cv","ridge_boot",
          "gam_k4","gam_k7","tree"),
  train=c(mean((tr$LC50 - p1)^2), mean((tr$LC50 - p2)^2),
          mean((tr$LC50 - predict(b_aic,tr))^2),
          mean((tr$LC50 - predict(f_aic,tr))^2),
          mean((tr$LC50 - predict(b_bic,tr))^2),
          mean((tr$LC50 - predict(f_bic,tr))^2),
          tr_cv, tr_boot, tr_g1, tr_g2, tr_tree),
  test=c(mean((te$LC50 - q1)^2), mean((te$LC50 - q2)^2),
         mean((te$LC50 - predict(b_aic,te))^2),
         mean((te$LC50 - predict(f_aic,te))^2),
         mean((te$LC50 - predict(b_bic,te))^2),
         mean((te$LC50 - predict(f_bic,te))^2),
         te_cv, te_boot, te_g1, te_g2, te_tree)
) %>% arrange(test)

print(res)
