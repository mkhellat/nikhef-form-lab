/*
 * Dump a fixed CDC-60 case vector.
 *
 * Reads validation/cases.txt. Each result line is
 *
 *   id op oct20 hex15 decimal
 *
 * under a "CDC_FP_ORACLE 1" header. The octal word is the low 60
 * bits, one's complement included, so a later FTN5 dump can be
 * compared without going through decimal conversion.
 *
 *   make -C cdc_fp oracle
 */
#include "cdc_fp.h"

#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_CASES 128
#define ID_LEN 8

struct item {
	char id[ID_LEN];
	char op[8];
	uint64_t w;
};

static void die(const char *msg)
{
	fprintf(stderr, "oracle_cdc_fp: %s\n", msg);
	exit(1);
}

static void die_line(int lineno, const char *msg)
{
	fprintf(stderr, "oracle_cdc_fp: line %d: %s\n", lineno, msg);
	exit(1);
}

static int split_tokens(char *line, char **tok, int max)
{
	int n = 0;
	char *p = line;

	while (*p) {
		while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n')
			p++;
		if (*p == '\0' || *p == '#')
			break;
		if (n == max)
			return -1;
		tok[n++] = p;
		while (*p && *p != ' ' && *p != '\t' && *p != '\r' &&
		       *p != '\n' && *p != '#')
			p++;
		if (*p == '#') {
			*p = '\0';
			break;
		}
		if (*p)
			*p++ = '\0';
	}
	return n;
}

static int valid_id(const char *id)
{
	size_t n = strlen(id);
	size_t i;

	if (n < 1 || n >= ID_LEN)
		return 0;
	if (!isalpha((unsigned char)id[0]))
		return 0;
	for (i = 1; i < n; i++) {
		if (!isalnum((unsigned char)id[i]))
			return 0;
	}
	return 1;
}

static int find_id(struct item *items, int n, const char *id, uint64_t *w)
{
	int i;

	for (i = 0; i < n; i++) {
		if (strcmp(items[i].id, id) == 0) {
			*w = items[i].w;
			return 0;
		}
	}
	return -1;
}

static int parse_lit(const char *tok, uint64_t *w)
{
	char *end;
	double x;

	if (strcmp(tok, "+0") == 0 || strcmp(tok, "0") == 0 ||
	    strcmp(tok, "+0.0") == 0) {
		*w = cdc_from_double(0.0);
		return 0;
	}
	if (strcmp(tok, "-0") == 0 || strcmp(tok, "-0.0") == 0) {
		*w = cdc_from_double(-0.0);
		return 0;
	}
	errno = 0;
	x = strtod(tok, &end);
	if (end == tok || *end != '\0' || errno == ERANGE)
		return -1;
	*w = cdc_from_double(x);
	return 0;
}

static void format_dec(uint64_t w, char *buf, size_t n)
{
	double d;

	if (cdc_is_indef(w)) {
		snprintf(buf, n, "%s",
			 (w & CDC_FP_SIGN) ? "-indef" : "+indef");
		return;
	}
	d = cdc_to_double(w);
	if (d == 0.0) {
		snprintf(buf, n, "%s", signbit(d) ? "-0" : "+0");
		return;
	}
	snprintf(buf, n, "%.16e", d);
}

static void emit(const struct item *it)
{
	char dec[64];

	format_dec(it->w, dec, sizeof dec);
	printf("%s %s %020" PRIo64 " %015" PRIx64 " %s\n", it->id, it->op,
	       it->w, it->w, dec);
}

int main(int argc, char **argv)
{
	FILE *fp;
	char line[256];
	struct item items[MAX_CASES];
	int n = 0;
	int lineno = 0;

	if (argc != 2)
		die("usage: oracle_cdc_fp cases.txt");
	fp = fopen(argv[1], "r");
	if (!fp) {
		fprintf(stderr, "oracle_cdc_fp: cannot open %s\n", argv[1]);
		return 1;
	}

	printf("CDC_FP_ORACLE 1\n");
	printf("# id op result_oct20 result_hex15 result_dec\n");

	while (fgets(line, sizeof line, fp)) {
		char *tok[5];
		int nt;
		struct item it;
		uint64_t a, b;
		long count, k;

		lineno++;
		if (strchr(line, '\n') == NULL && !feof(fp))
			die_line(lineno, "line longer than 255 characters");
		nt = split_tokens(line, tok, 5);
		if (nt == 0)
			continue;
		if (nt < 0)
			die_line(lineno, "too many tokens");
		if (!valid_id(tok[0]))
			die_line(lineno, "id must be 1..7 alphanumeric");
		memset(&it, 0, sizeof it);
		snprintf(it.id, sizeof it.id, "%s", tok[0]);
		if (find_id(items, n, it.id, &a) == 0)
			die_line(lineno, "duplicate id");

		if (nt == 3 && strcmp(tok[1], "lit") == 0) {
			snprintf(it.op, sizeof it.op, "lit");
			if (parse_lit(tok[2], &it.w) != 0)
				die_line(lineno, "bad lit");
		} else if (nt == 4 && (strcmp(tok[1], "add") == 0 ||
				       strcmp(tok[1], "sub") == 0 ||
				       strcmp(tok[1], "mul") == 0 ||
				       strcmp(tok[1], "div") == 0)) {
			snprintf(it.op, sizeof it.op, "%s", tok[1]);
			if (find_id(items, n, tok[2], &a) != 0 ||
			    find_id(items, n, tok[3], &b) != 0)
				die_line(lineno, "unknown operand id");
			if (strcmp(it.op, "add") == 0)
				it.w = cdc_add(a, b);
			else if (strcmp(it.op, "sub") == 0)
				it.w = cdc_sub(a, b);
			else if (strcmp(it.op, "mul") == 0)
				it.w = cdc_mul(a, b);
			else
				it.w = cdc_div(a, b);
		} else if (nt == 5 && strcmp(tok[1], "muln") == 0) {
			char *end;

			snprintf(it.op, sizeof it.op, "muln");
			if (find_id(items, n, tok[2], &a) != 0 ||
			    find_id(items, n, tok[3], &b) != 0)
				die_line(lineno, "unknown operand id");
			errno = 0;
			count = strtol(tok[4], &end, 10);
			if (end == tok[4] || *end != '\0' || errno != 0 ||
			    count < 1 || count > 64)
				die_line(lineno, "muln count must be 1..64");
			it.w = a;
			for (k = 0; k < count; k++)
				it.w = cdc_mul(it.w, b);
		} else {
			die_line(lineno, "want lit, add, sub, mul, div, or muln");
		}

		if (n >= MAX_CASES)
			die_line(lineno, "too many cases");
		items[n++] = it;
		emit(&it);
	}
	if (ferror(fp))
		die("read error");
	fclose(fp);
	if (n == 0)
		die("no cases");
	return 0;
}
