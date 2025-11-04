
library(tidyverse)
library(MASS)
library(rpart)
library(rpart.plot)
library(glmnet)
library(mgcv)
set.seed(4300)
theme_set(theme_bw())

# problem 1

df <- read.csv("data/qsar_aquatic_toxicity.csv",sep = ";", header = FALSE) %>% 
  as_tibble()

# a)

colnames(df) <- c("TPSA","SAacc","H050","MLOGP","RDCHI","GATS1p","nN","C040","LC50")

i <- sample(1:nrow(df), round(2 * nrow(df) / 3))
tr <- df[i, ]
te <- df[-i, ]

m1 <- lm(LC50 ~ ., data = tr)
m2 <- lm(LC50 ~ TPSA + SAacc + MLOGP + RDCHI + GATS1p + I(H050 > 0) + I(nN >
                                                                          0) + I(C040 > 0),
         data = tr)

p1 <- predict(m1, tr)
q1 <- predict(m1, te)
p2 <- predict(m2, tr)
q2 <- predict(m2, te)

cat("train_mse_m1", mean((tr$LC50 - p1)^2), "test_mse_m1", mean((te$LC50 -
                                                                   q1)^2), "\n")
cat("train_mse_m2", mean((tr$LC50 - p2)^2), "test_mse_m2", mean((te$LC50 -
                                                                   q2)^2), "\n")

p <- ggplot(tibble(obs = te$LC50, pred = q1), aes(obs, pred)) + geom_point() +
  geom_abline()
ggsave("plots/plot1a.pdf", p, width = 10, height = 6)


# b)

e1 <- c()
e2 <- c()

for (k in 1:200) {
  i <- sample(1:nrow(df), round(2 * nrow(df) / 3))
  tr <- df[i, ]
  te <- df[-i, ]
  m1 <- lm(LC50 ~ ., data = tr)
  m2 <- lm(LC50 ~ TPSA + SAacc + MLOGP + RDCHI + GATS1p + I(H050 > 0) +
             I(nN > 0) + I(C040 > 0),
           data = tr)
  e1[k] <- mean((te$LC50 - predict(m1, te))^2)
  e2[k] <- mean((te$LC50 - predict(m2, te))^2)
}

t <- tibble(model = c(rep("raw", 200), rep("dummy", 200)), mse = c(e1, e2))

p <- ggplot(t, aes(model, mse)) + geom_boxplot()
ggsave("plots/plot1b.pdf", p, width = 10, height = 6)

cat("mean_test_mse_raw", mean(e1), "\n")
cat("mean_test_mse_dummy", mean(e2), "\n")

# c)

full <- lm(LC50 ~ ., data=tr)
null <- lm(LC50 ~ 1, data=tr)

b_aic <- step(full, direction="backward", k=2, trace=0)
f_aic <- step(null, scope=formula(full), direction="forward", k=2, trace=0)
b_bic <- step(full, direction="backward", k=log(nrow(tr)), trace=0)
f_bic <- step(null, scope=formula(full), direction="forward", k=log(nrow(tr)), trace=0)

cat("backward_AIC:", paste(deparse(formula(b_aic), width.cutoff=500), collapse=" "), "\n")
cat("forward_AIC :", paste(deparse(formula(f_aic), width.cutoff=500), collapse=" "), "\n")
cat("backward_BIC:", paste(deparse(formula(b_bic), width.cutoff=500), collapse=" "), "\n")
cat("forward_BIC :", paste(deparse(formula(f_bic), width.cutoff=500), collapse=" "), "\n")

# d)

xtr <- as.matrix(dplyr::select(tr, -LC50)); ytr <- tr$LC50
xte <- as.matrix(dplyr::select(te, -LC50)); yte <- te$LC50
lam <- 10^seq(-4,4,length=100)
cv <- cv.glmnet(xtr,ytr,alpha=0,lambda=lam,standardize=TRUE)
lam_cv <- cv$lambda.min
tr_cv <- mean((ytr - as.numeric(predict(cv$glmnet.fit,s=lam_cv,newx=xtr)))^2)
te_cv <- mean((yte - as.numeric(predict(cv$glmnet.fit,s=lam_cv,newx=xte)))^2)
set.seed(4300)
B <- 200
oob_mse <- matrix(NA,B,length(lam))
for(b in 1:B){
  idx <- sample(seq_len(nrow(tr)), replace=TRUE)
  oob <- setdiff(seq_len(nrow(tr)), unique(idx))
  fit <- glmnet(xtr[idx,], ytr[idx], alpha=0, lambda=lam, standardize=TRUE)
  if(length(oob)>5){
    pr <- predict(fit, newx=xtr[oob,])
    oob_mse[b,] <- colMeans((matrix(ytr[oob],nrow=length(oob),ncol=length(lam)) - pr)^2)
  }
}
boot_curve <- colMeans(oob_mse, na.rm=TRUE)
lam_boot <- lam[which.min(boot_curve)]
fit_boot <- glmnet(xtr,ytr,alpha=0,lambda=lam,standardize=TRUE)
tr_boot <- mean((ytr - as.numeric(predict(fit_boot,s=lam_boot,newx=xtr)))^2)
te_boot <- mean((yte - as.numeric(predict(fit_boot,s=lam_boot,newx=xte)))^2)
p <- tibble(log_lambda=log(lam),
            cv=cv$cvm[match(lam, cv$lambda)],
            boot=boot_curve) %>% 
  pivot_longer(-log_lambda, names_to="proc", values_to="mse") %>% 
  ggplot(aes(log_lambda,mse,linetype=proc))+geom_line()

ggsave("plots/plot1d.pdf", p, width=10, height=6)

# e)

g1 <- gam(
  LC50 ~ s(TPSA,k=4) + s(SAacc,k=4) + s(H050,k=3) + s(MLOGP,k=4) +
    s(RDCHI,k=4) + s(GATS1p,k=4) + s(nN,k=3) + s(C040,k=3),
  data=tr, method="REML")

g2 <- gam(
  LC50 ~ s(TPSA,k=7) + s(SAacc,k=7) + s(H050,k=3) + s(MLOGP,k=7) +
    s(RDCHI,k=7) + s(GATS1p,k=7) + s(nN,k=3) + s(C040,k=3),
  data=tr, method="REML")

tr_g1 <- mean((tr$LC50 - predict(g1,tr))^2)
te_g1 <- mean((te$LC50 - predict(g1,te))^2)
tr_g2 <- mean((tr$LC50 - predict(g2,tr))^2)
te_g2 <- mean((te$LC50 - predict(g2,te))^2)

# f)

tree0 <- rpart(LC50~., data=tr, method="anova", cp=0.001)
cp_opt <- tree0$cptable[which.min(tree0$cptable[,"xerror"]), "CP"]
tree <- prune(tree0, cp=cp_opt)
tr_tree <- mean((tr$LC50 - predict(tree,tr))^2); te_tree <- mean((te$LC50 - predict(tree,te))^2)
pdf("plots/plot1f_tree.pdf", width=10, height=6); rpart.plot(tree); dev.off()
pdf("plots/plot1f_cp.pdf", width=10, height=6); plotcp(tree0); dev.off()

# g)

MSE_tr_m1 <- mean((tr$LC50 - p1)^2)
MSE_te_m1 <- mean((te$LC50 - q1)^2)
MSE_tr_m2 <- mean((tr$LC50 - p2)^2)
MSE_te_m2 <- mean((te$LC50 - q2)^2)
tr_ba <- mean((tr$LC50 - predict(b_aic,tr))^2); te_ba <- mean((te$LC50 - predict(b_aic,te))^2)
tr_fa <- mean((tr$LC50 - predict(f_aic,tr))^2); te_fa <- mean((te$LC50 - predict(f_aic,te))^2)
tr_bb <- mean((tr$LC50 - predict(b_bic,tr))^2); te_bb <- mean((te$LC50 - predict(b_bic,te))^2)
tr_fb <- mean((tr$LC50 - predict(f_bic,tr))^2); te_fb <- mean((te$LC50 - predict(f_bic,te))^2)

res <- tibble(
  model=c("lm_raw","lm_dummy","step_back_AIC","step_fwd_AIC","step_back_BIC","step_fwd_BIC","ridge_cv","ridge_boot","gam_k4","gam_k8","tree"),
  train=c(MSE_tr_m1,MSE_tr_m2,tr_ba,tr_fa,tr_bb,tr_fb,tr_cv,tr_boot,tr_g1,tr_g2,tr_tree),
  test =c(MSE_te_m1,MSE_te_m2,te_ba,te_fa,te_bb,te_fb,te_cv,te_boot,te_g1,te_g2,te_tree)
) %>% arrange(test)
print(res)

