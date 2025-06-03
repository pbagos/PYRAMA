import meta_analysis
import beta_se_meta2
import sys

def robust_methods(data, output_folder, base_filename, inheritance_model, effect_size_type,
                   robust_method, type_of_effect, approximate_max, biv_ma='NO'):
    print("Applying Robust Methods:")
    case_1_columns = ['SNP', 'CHR', 'BP', 'aa1', 'ab1', 'bb1', 'aa0', 'ab0', 'bb0']
    case_2_columns = ['SNP', 'CHR', 'BP', 'BETA', 'SE']

    if all(col in data.columns for col in case_1_columns):
        data_subset = data[case_1_columns]
        result = meta_analysis.meta_analysis(
            data_subset,
            inheritance_model='ALL',
            effect_size_type=effect_size_type,
            robust_method=robust_method,
            type_of_effect=type_of_effect,
            approximate_max=approximate_max
        )
    elif all(col in data.columns for col in case_2_columns):
        data_subset = data[case_2_columns]
        result = beta_se_meta2.beta_se(data_subset, biv_beta_input=biv_ma)
    else:
        raise ValueError("Data does not match the required columns for Case 1 or Case 2.")

    result.to_csv(f"{output_folder}/{base_filename}_robust_methods_results.txt", sep='\t', index=False)
    print(f"Robust methods results saved to {output_folder}/{base_filename}_robust_methods_results.txt")
    return result


if __name__ == "__main__":
    robust_methods(*sys.argv[1:])

