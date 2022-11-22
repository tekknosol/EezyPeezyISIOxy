# EezyPeezyISIOxy

This repository is being used to develop a workflow for analyzing the ISIMIP-3a data to study the effects of climate change on lake hypolimnetic oxygen dynamics and originated at the GLEON 2022 conference at Lake George, NY.

Contacts: co-champions Carolina Barbosa and Philipp Keller

Working group members: Cayelan Carey, Chloe Faehndrich, Robert Ladwig, Sofia LaFuente, Rafa Marce, Daniel Mercado-Battin, Lipa Nwala

GIT-HUB sources:
https://www.youtube.com/watch?v=B-FHx4l1BNU 

<a href="url"><img src="OxygenTest.png" width=80% height=80% ></a>

## Flowchart

```mermaid
graph TD
    A[1 Folder per Lake - Lake 1:n] --> B(output_temp.txt)
    A --> C(output_z.txt)
    A --> D(hypsograph.dat)
    subgraph ide1 [" "]
    B & C & D
    end
    B --> E{Thermal script}
    C --> E
    D --> E
    E --> F(thermal_info.csv)
    F --> G{Oxygen script}
    G --> H(oxygen_info.csv)
    F --to be decided--> I(output.nc)
    H --to be decided--> I

```
