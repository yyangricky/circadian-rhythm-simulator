# circadian-rhythm-simulator
Simulates the effects of sleep loss according to established circadian rhythm models.

This code uses the `circadian` package. You do not need to install it, although it can be pip stalled using `pip install circadian`. The latest python
version is recommended.

## How to Run

### Running `graphactogram.py`
This script generates actograms based on circadian rhythm models. To run it:
1. Ensure you have Python installed.
2. Execute the script:
    ```bash
    python graphactogram.py
    ```
3. The script will display the generated actogram. For automated runs, consider saving the plot:
    ```bash
    python graphactogram.py --save
    ```

### Running `improvedfigure.py`
This script reproduces enhanced figures based on the circadian rhythm simulations. To run it:
1. Ensure all dependencies are installed.
2. Execute the script:
    ```bash
    python improvedfigure.py
    ```
3. The output will be displayed or saved depending on the script's configuration.

For both scripts, you can modify parameters directly in the code or pass arguments if supported. Refer to the script comments for customization options.