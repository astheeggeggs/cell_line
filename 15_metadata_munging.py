import hail as hl
hl.init()

IMPUTESEX_FILE = '/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/celldataB37_dp_impute_sex.tsv'
Y_NCALLED = '/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/celldataB37_dp_y_called.tsv'

# Read in plink file and run sex imputation
X_PRUNED = '/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/celldataB37_dp_X_pruned'
PLINK_CELL = '/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/celldataB37_dp'

mt = hl.import_plink(
	bed = X_PRUNED + '.bed',
	bim = X_PRUNED + '.bim',
	fam = X_PRUNED + '.fam'
	)

mt_X = mt.filter_rows(mt.locus.in_x_nonpar())
imputed_sex = hl.impute_sex(mt_X.GT)
mt = mt.annotate_cols(impute_sex = imputed_sex[mt.s])

mt.cols().select('impute_sex').flatten().export(IMPUTESEX_FILE)

# Determine non-missing allele count on the y.
mt = hl.import_plink(
	bed = PLINK_CELL + '.bed',
	bim = PLINK_CELL + '.bim',
	fam = PLINK_CELL + '.fam'
	)

mt = mt.filter_rows(mt.locus.in_y_nonpar() | mt.locus.in_y_par())
mt = hl.sample_qc(mt, name='qc')

mt = mt.cols()
mt.select(n_called=mt.qc.n_called).export(Y_NCALLED)

# Also, evaluate the relatedness of the sequences in the cell-lines. Are the replicates really 100% related?
PLINK_CELL = '/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/celldataB37_dp_test-updated-autosomes-combined_pruned'
mt = hl.import_plink(
	bed = PLINK_CELL + '.bed',
	bim = PLINK_CELL + '.bim',
	fam = PLINK_CELL + '.fam'
	)

IBD_OUTPUT = "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/celldataB37_ibd.tsv"
ibd_table = hl.identity_by_descent(mt)
ibd_table = ibd_table.flatten()
ibd_table.export(IBD_OUTPUT)
