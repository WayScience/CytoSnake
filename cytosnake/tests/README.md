# CytoSnake Testing

- [CytoSnake Testing](#cytosnake-testing)
  - [Rationale and Motivation](#rationale-and-motivation)
    - [Performance](#performance)
    - [Reproducibility](#reproducibility)
    - [User-friendliness and Portability](#user-friendliness-and-portability)
  - [Best Practices](#best-practices)
    - [Tech stack](#tech-stack)
    - [Dataset](#dataset)
    - [Clear and Descriptive test](#clear-and-descriptive-test)
  - [Testing Module Main Docuemtnation.](#testing-module-main-docuemtnation)

---

Welcome to the testing documentation for `CytosSnake`, where you'll learn the essential principles of crafting effective and insightful tests.
Our tests focus on

In this documentation, you will learn the mission of making effective tests in `CytoSnake` Ensuring that every user who uses `CytoSnake` has

## Rationale and Motivation

Below are the rationals and motivation behind generating tests for `CytoSnake`

### Performance 

Remarkable technological advances have occurred within the field of high-throughput microscopy imaging, generating a surmount of single-cell morphology readouts, also known as single-cell morphology profiles. 
With this increase in data, there is a high need for high-performing software tools that can not only scale with large volumes of single-cell morphology profiles but also analyze the quickly.  
Our tests in `CytoSnake` are designed to align with our goals of scale and performance. By rigorously evaluating core functionality, we ensure that users can quickly obtain conclusions from their data, enabling faster insights and analysis.

### Reproducibility

Reproducibility is an integral part of Cytosnake ensuring the reliable robust conclusions 

Our tests are designed to ensure reproducibility by monitoring `CytoSnake`'s functionality across different scenarios.
This means that `CytoSnake`'s will be subjected to various inputs and conditions to test actively identifying and capturing any inconsistencies with the generated outputs. 
Our testing design aims to offer researchers in the community a sense of reliability, ensuring they can consistently achieve reproducible results that reinforce robust scientific conclusions.
Overall, `CytoSnake`'s tests are designed to not only reach robust conclusions but also accelerate the progress in the field of Cell biology, enabling faster and more reliable advancements. 


### User-friendliness and Portability

While prioritizing reproducibility is crucial, we also intend to make `CytoSnake` accessible to a broader audience, ensuring its user-friendliness and portability.
Our tests are deployed in GitHub Actions, allowing for multi-platform testing acorss different operating systems, ensuring `CytoSnake` is system agnostic.
In addition, we use `Poetry` package to manage dependencies, assuring that the software can be easily installed.
We leverage the `conda` environment manager to encapsulate dependencies within `CytoSnake`'s modules allowing it to be portable and consistent across all systems. 
Ultimately, our tests aim to provide users with an easy and enjoyable experience while using `CytoSnake`, by providing the convenience of being portable and requiring little effort to use it.

## Best Practices

We believe that implementing good testing practices can create software that's not only robust but also significantly more dependable.
This ensures a high level of reliability for users and developers to conduct their analysis with confidence in `CytoSnake`.

Outlined below are the best practices that we consider instrumental in helping us achieve our goals.

### Tech stack

Our testing suite uses [`pytest`] a very powerful and flexible testing

We also use `Code Coverage` as a metric for our tests, 

### Dataset

The datasets that we use are structured, emulating

### Clear and Descriptive test

In `CytoSnake`, well-structured testing documentation plays an important role.
It simplifies the software experience for both users and developers by providing a clear understanding of what is being tested within `CytoSnake`.
We believe that providing easy-to-follow documentation generated complete t

## Testing Module Main Docuemtnation. 

Below are the documetations for each testing suite:

- [Functional testing](./functional/README.md)
- [Unit testing]()
- [Workflow testing]()
