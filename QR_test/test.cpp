#include "pch.h"
#include "../QR-decomposition/QR_header.h"

//TEST(QR, testing_tests) {
//  EXPECT_EQ(1, 1);
//  EXPECT_TRUE(true);
//}

TEST(QR, is_correct_128_size_with_block_32) {
	QR <double> t(0, 128, 32);
	t.HHolder_A();
	t.HHolder_R();
	t.HHolder_Q();
	t.transpQ();
	EXPECT_TRUE(t.check());
}

TEST(QR, is_correct_512_size_with_block_48) {
	QR <double> t(0, 512, 48);
	t.HHolder_A();
	t.HHolder_R();
	t.HHolder_Q();
	t.transpQ();
	EXPECT_TRUE(t.check());
}

TEST(QR, is_correct_512_size_with_block_64) {
	QR <double> t(0, 512, 64);
	t.HHolder_A();
	t.HHolder_R();
	t.HHolder_Q();
	t.transpQ();
	EXPECT_TRUE(t.check());
}

TEST(QR, is_correct_1024_size_with_block_64) {
	QR <double> t(0, 1024, 64);
	t.HHolder_A();
	t.HHolder_R();
	t.HHolder_Q();
	t.transpQ();
	EXPECT_TRUE(t.check());
}