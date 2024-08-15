# SiteKIT

SiteKIT can rapidly and accurately identify key sites that may be linked to specific traits. Site-trait association analysis aim to identify associations of biological traits with sites mutations. These mutations may affect gene function, thereby altering the phenotypic characteristics of organisms or causing diseases. SiteKIT provides an opportunity to gain insights into the complex relationship between site mutations and biological traits.

# Installation

There are two ways to use SiteKIT:

1. **Using the standalone executable (SiteKIT):**
   - This approach involves using a pre-packaged executable that contains all necessary dependencies, making it easy to run without needing to set up a specific environment.
2. **Configuring environment dependencies via `requirements.txt` and using the Python script (SiteKIT.py):**
   - In this method, you configure the required packages by listing them in a `requirements.txt` file. Then, you can run the Python script (SiteKIT.py) within an environment that has these dependencies installed. The version of Python is 3.12.5

# Tutorial

### Parameters

The `SiteKIT` tool supports various parameters to customize its behavior. You can view the help instructions by running `SiteKIT -h`.

- **`--aa_dir`**: The required input folder containing multiple sequence alignments (MSAs) of amino acids.
- **`--group_id`**: Specifies IDs or prefixes for trait groups. Include the trait group prefix list in single quotes ('') separated by commas (,). Do not include spaces. For example: `'sample1,sample2,...'`.
- **`--codon_dir`**: The input folder containing multiple sequence alignments (MSAs) of codons.
- **`-o`, `--output_prefix`**: The name of the output folder. [Default: `Site_Trait_Association_Analysis`]
- `-m`, `--mode`:
  - **`mode1`**: Preserves the sequence ID exactly as it is. [Default]
  - **`mode2`**: Replaces any non-alphanumeric characters in the sequence ID with '_' (underscore), while retaining the content before the first space.
- **`--separator`**: Not enabled by default. Specify the separator when using `mode2` to enable it.
- **`--multiprocessing`**: The number of processes to use for parallel processing. [Default: 10]

### Example

Within the `Example_Marsupials` directory, you can execute the following commands to test the AA (amino acid) and codon modes:


```commandline
#AA mode
SiteKIT --aa_dir AA  --group_id Drom,Phas,Macr -o Marsupials_AA --mode mode2 --separator ':'
#codon mode
SiteKIT --aa_dir AA  --codon_dir Codon --group_id Drom,Phas,Macr -o Marsupials_Codon --mode mode2 --separator ':'
```

In `Example_Penguins`, the following command can be executed for testing.

```
SiteKIT --aa_dir Sphenisciformes_genes  --group_id Sphenisciformes -o Sphenisciformes_AA
```



# Citations



