import unittest
from pathlib import Path

from msspack.config import (
    InputsConfig,
    MSSPackConfig,
    PipelineConfig,
    ProjectConfig,
    ReferenceConfig,
    SampleConfig,
    StCommentConfig,
    SubmissionConfig,
    SubmitterConfig,
    ToolsConfig,
)
from msspack.header import render_header
from msspack.submission_render import render_final_annotation


def _config() -> MSSPackConfig:
    return MSSPackConfig(
        base_dir=Path("/tmp"),
        project=ProjectConfig(name="Demo"),
        inputs=InputsConfig(fasta="a.fa", gff="a.gff"),
        sample=SampleConfig(
            locus_tag="Demo",
            locus_tag_digits=6,
            scientific_name="Demo species",
            isolate="X1",
            geo_loc_name="Japan",
            collection_date="2026-01-01",
        ),
        submission=SubmissionConfig(
            datatype="WGS",
            hold_date="20261231",
            bioproject="PRJDB1",
            biosample="SAMD1",
            sra=["DRR1", "DRR2"],
            keywords=["WGS", "STANDARD_DRAFT"],
        ),
        submitter=SubmitterConfig(
            ab_name=["Fukushima,K."],
            contact="Kenji Fukushima",
            institute="NIG",
            department="Lab",
            country="Japan",
            state="Shizuoka",
            city="Mishima",
            street="1111 Yata",
            zip="411-8540",
            phone="81-00-0000-0000",
            email="x@example.org",
        ),
        reference=ReferenceConfig(
            title="Demo sequencing",
            ab_name=["Fukushima,K."],
            year=2026,
        ),
        st_comment=StCommentConfig(
            assembly_method="Assembler 1.0",
            assembly_name="Demo_1.0",
            genome_coverage="100X",
            sequencing_technology="Nanopore",
        ),
        pipeline=PipelineConfig(),
        tools=ToolsConfig(),
    )


class HeaderRenderTests(unittest.TestCase):
    def test_render_header_contains_core_sections(self) -> None:
        text = render_header(_config())
        self.assertIn("COMMON\tDATE\t\thold_date\t20261231", text)
        self.assertIn("\tDBLINK\t\tproject\tPRJDB1", text)
        self.assertIn("\t\t\tbiosample\tSAMD1", text)
        self.assertIn("\tKEYWORD\t\tkeyword\tWGS", text)
        self.assertIn("\tSUBMITTER\t\tab_name\tFukushima,K.", text)
        self.assertIn("\tREFERENCE\t\ttitle\tDemo sequencing", text)
        self.assertIn("\tST_COMMENT\t\ttagset_id\tGenome-Assembly-Data", text)

    def test_render_header_omits_date_for_immediate_release(self) -> None:
        config = _config()
        config.submission.hold_date = ""

        text = render_header(config)

        self.assertTrue(text.startswith("COMMON\tDATATYPE\t\ttype\tWGS\n"))
        self.assertNotIn("\tDATE\t", text)
        self.assertNotIn("hold_date", text)

    def test_render_final_annotation_normalizes_country_to_geo_loc_name(self) -> None:
        final = render_final_annotation(
            "COMMON\tDATE\t\thold_date\t20261231\n",
            "\tCDS\t1..3\tlocus_tag\tDemo_000001\n\t\t\tcountry\tJapan\n",
        )
        self.assertIn("\t\t\tgeo_loc_name\tJapan\n", final)
        self.assertNotIn("\t\t\tcountry\tJapan\n", final)
