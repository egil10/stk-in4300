# STK-IN4300 – Statistical Learning Methods in Data Science

Repository containing assignments and resources for **STK-IN4300 – Statistical Learning Methods in Data Science** at the University of Oslo.

**Course Link**: [STK-IN4300 Course Page](https://www.uio.no/studier/emner/matnat/math/STK-IN4300/index-eng.html)

## Course Description

This course focuses on methods for modern data analysis both within a practical and theoretical framework. Fewer assumptions are made when using these methods, such as machine learning or statistical learning, than when using classical methods. Consequently, the data play a larger role in the identification of structures and relationships. Starting from the basic methods, the course will then cover more advanced procedures, specifically designed to tackle modern data challenges such as increasing complexity and large amounts of information (Big Data settings).

## Learning Outcomes

After completing the course you will:

- Understand key concepts for a good analysis of the data
- Understand the theoretical aspects of methods within machine/statistical learning
- Know a range of different methods for data analysis, including:
  - Penalized likelihood and basis expansions
  - Neural networks
  - Boosting and ensemble methods
  - Gaussian processes within machine learning
- Know procedures for fitting such methods to data, including:
  - (Stochastic) gradient descent
  - Back-propagation
- Be able to evaluate strengths and weaknesses of the methods and choose between them in practice

## Repository Structure

```
stk-in4300/
│
├── oblig1/                    # Assignment 1
│   ├── data/                  # Data files for assignment 1
│   ├── plots/                 # Generated plots and visualizations
│   │   ├── child_height_bad.pdf
│   │   ├── child_height_good.pdf
│   │   ├── parent_sex_distribution.pdf
│   │   └── parent_sex_distribution_good.pdf
│   ├── problem1.R            # R script for problem 1
│   ├── problem2.R            # R script for problem 2
│   └── STK-IN4300 Oblig 1.pdf # Assignment description
│
├── oblig2/                    # Assignment 2
│   ├── data/                  # Data files for assignment 2
│   │   └── qsar_aquatic_toxicity.csv
│   ├── plots/                 # Generated plots and visualizations
│   │   ├── plot1a.pdf
│   │   ├── plot1b.pdf
│   │   ├── plot1d.pdf
│   │   ├── plot1f_cp.pdf
│   │   ├── plot1f_tree.pdf
│   │   ├── plot2a.pdf
│   │   └── plot2c_tree.pdf
│   ├── 1.R                   # R script for problem 1
│   ├── 2.R                   # R script for problem 2
│   └── STK-IN4300 Oblig 2.pdf # Assignment description
│
├── rsc/                       # Resources
│   ├── stk-in4300 2023-2018.pdf       # Past exam papers (2018-2023)
│   ├── stk-in4300 oblig 1.pdf         # Assignment 1 resources
│   └── stk-in4300 oblig 2.pdf         # Assignment 2 resources
│
└── README.md                  # This file
```

## Assignments

### Assignment 1

Explores foundational statistical learning concepts, including:
- Data exploration and visualization
- Statistical modeling and analysis
- Use of the `HistData` package (PearsonLee dataset)

**Files:**
- `problem1.R` - Analysis scripts
- `problem2.R` - Additional problem solutions
- Generated plots in `plots/` directory

### Assignment 2

Focuses on advanced regression and classification methods:
- QSAR (Quantitative Structure-Activity Relationship) aquatic toxicity regression
- Tree-based methods (using `rpart`)
- Regularized regression methods (using `glmnet`)
- Generalized additive models (using `mgcv`)
- Cross-validation and model selection

**Files:**
- `1.R` - Problem 1 solutions
- `2.R` - Problem 2 solutions
- `data/qsar_aquatic_toxicity.csv` - Dataset for analysis

## Technologies Used

- **R** - Primary programming language
- **R Packages:**
  - `tidyverse` - Data manipulation and visualization
  - `MASS` - Statistical functions
  - `rpart` / `rpart.plot` - Decision trees
  - `glmnet` - Regularized regression
  - `mgcv` - Generalized additive models
  - `HistData` - Historical datasets

## Course Information

- **Credits**: 10
- **Level**: Master
- **Teaching**: Autumn
- **Examination**: Autumn
- **Language**: English

### Recommended Previous Knowledge

- Basic mathematical knowledge in probability (equal to STK1100)
- Linear regression (equal to STK1110)
- Linear algebra (equal to MAT1120 or MAT1125)
- Basic programming skills (equal to IN1900)

### Examination

- Final written or oral exam (100% of final grade)
- 2 mandatory assignments that must be approved before sitting the final exam
- No examination support material allowed

## Resources

- [Official Course Page](https://www.uio.no/studier/emner/matnat/math/STK-IN4300/index-eng.html)
- Past exam papers (2018-2023) available in `rsc/` directory

## Contact

For questions about the course, contact the **Department of Mathematics** at the University of Oslo.

---

*Last updated: December 2025*
