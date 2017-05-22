---
title: "SDSU"
layout: default
git_repo: "sdsu"
---
 
{% assign repo = site.data.git_repo[page.git_repo] %}

# South Dakota State University

## Git Repository
Found at [{{ repo.github_url }}]({{ repo.github_url }})

{% include_relative README.md %}

## Classes

{% for class in site.data.sdsu.classes %}
- [{{ class.uid }}](/sdsu/{{ class.uid }})
    - {{class.title}}
{% endfor %}
