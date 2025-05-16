# galaxytools


Collection of Galaxy tool wrappers maintained and developed by [ELIXIR Belgium](https://www.elixir-belgium.org/).


## General guide on coverting a Python tool into a Galaxy wrapper

WIP: example with protein_calculator:

- Create a virtual env for your tool, containing Planemo

```sh
python3 -m venv .
. bin/activate
```

- Install planemo and initiate the tool wrapper xlm file
```sh
pip install planemo
planemo tool_init --id 'protein_calculator' --name 'Protein Calculator from the VIB Protein Core'
``` 

- Initialise tool wrapper
```sh
planemo tool_init --force \
                    --id 'protein_calculator' \
                    --name 'Protein Calculator from the VIB Protein Core' \
                    --requirement python@3.9 \
                    --example_command 'python main.py --name 'test_prot1' --sequence "LLLLLLEEEEEVVVVV"' \
                    --example_input 'sequence.txt' \
                    --example_output 'report.html' \
                    --test_case \
                    --cite_url 'https://github.com/KBalcaen/Protein_Calculator' \
                    --help_from_command 'python main.py --help' \
		    --autopygen main.py
```
